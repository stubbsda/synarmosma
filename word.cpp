#include "word.h"

using namespace SYNARMOSMA;

extern Random RND;

Word::Word()
{

}

Word::Word(unsigned int n)
{
  initialize(n);
}

Word::Word(unsigned int n,int m)
{
  initialize(n,m);
}

Word::Word(const std::string& w)
{
  unsigned int i,j,n;
  int e;
  char alphabet[] = {'a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z'};
  Word temp;

  for(i=0; i<w.length(); ++i) {
    assert(std::isalpha(w[i]));
    for(j=0; j<26; ++j) {
      if (std::tolower(w[i]) == alphabet[j]) {
        n = j;
        break;
      }
    }
    e = (std::isupper(w[i])) ? -1 : 1;
    temp.content.push_back(std::pair<unsigned int,int>(n,e));
  }
  // This word needs to be normalized now...
  *this = temp.normalize();
}

Word::Word(const Word& source)
{
  content = source.content;
}

Word& Word::operator =(const Word& source)
{
  if (this == &source) return *this;

  content = source.content;

  return *this;
}

Word::~Word()
{

}

void Word::initialize(unsigned int n)
{
  content.clear();
  if (n == 0) return;

  unsigned int i;
  std::pair<unsigned int,int> doublet;

  for(i=0; i<n; ++i) {
    doublet.first = i; 
    doublet.second = RND.irandom(1,10);
    if (RND.drandom() < 0.5) doublet.second *= -1;
    content.push_back(doublet);
  }
}

void Word::initialize(unsigned int n,int m)
{
  if (m == 0) throw std::invalid_argument("The exponent of the word's only letter must not be zero!");

  std::pair<unsigned int,int> doublet(n,m);
  content.clear();
  content.push_back(doublet);
}

void Word::initialize(const std::vector<unsigned int>& base,const std::vector<int>& exponent)
{
  if (base.size() != exponent.size()) throw std::invalid_argument("The length of the base and exponent vectors for word initialization must be identical!");

  unsigned int i;
  std::pair<unsigned int,int> doublet;

  content.clear();
  for(i=0; i<base.size(); ++i) {
    if (exponent[i] == 0) continue;
    doublet.first = base[i];
    doublet.second = exponent[i];
    content.push_back(doublet);
  }
}

unsigned int Word::get_alphabet(std::set<unsigned int>& alphabet) const
{
  // This method computes the distinct letters in this word...
  std::vector<std::pair<unsigned int,int> >::const_iterator it;

  alphabet.clear();
  for(it=content.begin(); it!=content.end(); ++it) {
    alphabet.insert(it->first);
  }
  return alphabet.size();
}

bool Word::legal() const
{
  std::vector<std::pair<unsigned int,int> >::const_iterator it;

  for(it=content.begin(); it!=content.end(); ++it) {
    if (it->second == 0) return false;
  }
  return true;
}

bool Word::trivial() const
{
  if (content.size() == 1) {
    if (content[0].second == 1 || content[0].second == -1) return true;
  }
  return false;
}

bool Word::alias() const
{
  if (content.size() == 2) {
    if (content[0].first == content[1].first) throw std::runtime_error("The class instance in Word::alias has not been normalized!");
    if ((std::abs(content[0].second) == 1) && (std::abs(content[1].second) == 1)) return true;
  }
  return false;
}

bool Word::homogeneous() const
{
  // This is needed for the algorithm to test words in a braid group for their 
  // equivalence
  unsigned int prime;
  bool ppow = false,npow = false;
  std::set<unsigned int> S;
  std::vector<std::pair<unsigned int,int> >::const_iterator it;

  // Get the smallest letter in this word's alphabet...
  get_alphabet(S);
  prime = *S.begin();

  for(it=content.begin(); it!=content.end(); ++it) {
    if (it->first == prime) {
      if (it->second > 0) {
        ppow = true;
      }
      else if (it->second < 0) {
        npow = true;
      }
    }
  }
  if (ppow && npow) return false;
  return true;
}

Word Word::permute(unsigned int n) const
{
  // This method will create the word formed by the cyclic
  // permutation of the current instance, that is
  // a_1*a_2*...*a_k -> a_n*a_{n+1}*...*a_k*a_1*a_2*...*a_{n-1}
  if (n >= length()) throw std::invalid_argument("The permutation index must be less than the word's length!");
  Word output;

  if (n == 0) {
    output = Word(*this);
    return output;
  }
  unsigned int i;
  for(i=n; i<length(); ++i) {
    output.content.push_back(content[i]);
  }
  for(i=0; i<n; ++i) {
    output.content.push_back(content[i]);
  }
  return output;
}

Word Word::operator !() const
{
  return (this->invert());
}

Word Word::invert() const
{  
  int i,n = (signed) content.size();
  Word output; 
  std::pair<unsigned int,int> doublet;

  for(i=n-1; i>=0; --i) {
    doublet.first = content[i].first;
    doublet.second = -content[i].second;
    output.content.push_back(doublet);
  }
  return output;
}

Word Word::swap(unsigned int p,unsigned int q,bool invert) const
{
  unsigned int i,n = content.size();
  Word output;
  std::pair<unsigned int,int> doublet;

  for(i=0; i<n; ++i) {
    doublet.first = content[i].first;
    doublet.second = content[i].second;
    if (doublet.first == q) {
      doublet.first = p;
      if (invert) doublet.second *= -1;
    }
    else if (doublet.first > q) {
      doublet.first -= 1;
    }
    output.content.push_back(doublet);
  }
  return output;
}

Word Word::mutate() const
{
  if (empty()) throw std::runtime_error("An empty word cannot be mutated!");

  Word output;
  std::set<unsigned int> S;
  std::pair<unsigned int,int> doublet;
  int n = RND.irandom(length()); 
  unsigned int m = get_alphabet(S);

  output.content = content;

  if (m == 1) {
#ifdef DEBUG
    assert(n == 1);
#endif
    output.content[0].first = content[0].first;
    do {
      n =  RND.irandom(1,10);
      if (RND.drandom() < 0.5) n *= 1;
      if (n != content[0].second) break;
    } while(true);
    output.content[0].second = n;
  }
  else {
    do {
      doublet.first = RND.irandom(S); 
      if (doublet.first != content[n].first) break;
    } while(true);
    doublet.second = RND.irandom(1,10);
    if (RND.drandom() < 0.5) doublet.second = -doublet.second;
    output.content[n] = doublet;
  }

  return output;
}

int Word::serialize(std::ofstream& s) const
{
  unsigned int i,j,n = content.size();
  int k,count = 0;

  s.write((char*)(&n),sizeof(int)); count += sizeof(int);
  for(i=0; i<n; ++i) {
    j = content[i].first;
    s.write((char*)(&j),sizeof(int)); count += sizeof(int);
    k = content[i].second;
    s.write((char*)(&k),sizeof(int)); count += sizeof(int);
  }
  return count;
}

int Word::deserialize(std::ifstream& s)
{
  unsigned int i,n,base;
  int exponent,count = 0;
  std::pair<unsigned int,int> doublet;

  clear();

  s.read((char*)(&n),sizeof(int)); count += sizeof(int);
  for(i=0; i<n; ++i) {
    s.read((char*)(&base),sizeof(int)); count += sizeof(int);
    s.read((char*)(&exponent),sizeof(int)); count += sizeof(int);
    doublet.first = base;
    doublet.second = exponent;
    content.push_back(doublet);
  }
  return count;
}

Word Word::reduce() const
{
  unsigned int i,n,m,rpoint,L;
  bool altered;
  std::vector<std::pair<unsigned int,int> > vx;
  Word output;

  output.content = content;
  do {
    altered = false;
    L = output.length();
    for(i=0; i<L-1; ++i) {
      n = output.content[i].first;
      m = output.content[i+1].first;
      if (n == m) {
        if (output.content[i].second == -output.content[i+1].second) {
          rpoint = i;
          altered = true;
          break;
        }
      }
    }
    if (!altered) break;
    vx.clear();
    for(i=0; i<L; ++i) {
      if (i == rpoint || i == (rpoint+1)) continue;
      vx.push_back(output.content[i]);
    }
    output.content = vx;
  } while(true);

  return output;
}

Word Word::reduce(const std::set<unsigned int>& trivial_generators,const unsigned int* offset) const
{
  Word output;
  std::pair<unsigned int,int> doublet;
  std::vector<std::pair<unsigned int,int> >::const_iterator it;

  for(it=content.begin(); it!=content.end(); ++it) {
    if (trivial_generators.count(it->first) == 1) continue;
    doublet.first = offset[it->first];
    doublet.second = it->second;
    output.content.push_back(doublet);
  }
#ifdef DEBUG
  assert(output.legal());
#endif
  return output;
}

Word Word::normalize() const
{
  unsigned int i = 0;
  const unsigned int n = length();
  std::pair<unsigned int,int> doublet;
  Word output;

  if (n == 0) return output;

  if (n == 1) {
    if (content[0].second != 0) output.content = content;
    return output;
  }
  do {
    if (i >= n) break;
    if (i == (n - 1)) {
      if (content[i].second != 0) output.content.push_back(content[i]);
      break;
    }
    if (content[i].first == content[i+1].first) {
      doublet.first = content[i].first;
      doublet.second = content[i].second + content[i+1].second;
      if (doublet.second != 0) output.content.push_back(doublet);
      i += 2;
    }
    else {
      if (content[i].second != 0) output.content.push_back(content[i]);
      i += 1;
    }
  } while(true);

#ifdef DEBUG
  assert(output.legal());
#endif
  return output;
}

namespace SYNARMOSMA {
  bool operator ==(const Word& w1,const Word& w2)
  {
    unsigned int i;

    if (w1.content.size() == w2.content.size()) {
      bool eq = true;
      for(i=0; i<w1.content.size(); ++i) {
        if (w1.content[i] != w2.content[i]) {
          eq = false;
          break;
        }
      }
      if (eq) return true;
    }
    return false;
  }

  bool operator !=(const Word& w1,const Word& w2)
  {
    bool output = (w1 == w2) ? false : true;
    return output;
  }

  Word operator *(const Word& w1,const Word& w2)
  {
    unsigned int i;
    Word temp = w1;
 
    for(i=0; i<w2.content.size(); ++i) {
      temp.content.push_back(w2.content[i]);
    }
    Word output = temp.normalize();
    return output;
  }

  std::ostream& operator <<(std::ostream& os,const Word& source)
  {
    unsigned int i;
    std::pair<unsigned int,int> doublet;

    if (source.empty()) return os;

    for(i=0; i<source.content.size()-1; ++i) {
      doublet = source.content[i];
      if (doublet.second == 0) continue;
      if (doublet.second == 1) {
        os << "x[" << 1 + doublet.first << "]*";
      }
      else {
        os << "x[" << 1 + doublet.first << "]^(" << doublet.second << ")*";
      }
    }
    doublet = source.content[source.content.size()-1];
    if (doublet.second == 1) {
      os << "x[" << 1 + doublet.first << "]";
    }
    else {
      os << "x[" << 1 + doublet.first << "]^(" << doublet.second << ")";
    }
    return os;
  }

  int affinity(const Word& w1,const Word& w2,std::set<unsigned int>& output)
  {
    // This method computes the intersection of the alphabet of the two words
    unsigned int n;
    std::set<unsigned int> S1,S2;
    std::set<unsigned int>::const_iterator it;

    w1.get_alphabet(S1);
    w2.get_alphabet(S2);

    output.clear();

    for(it=S1.begin(); it!=S1.end(); ++it) {
      n = *it;
      if (S2.count(n) > 0) output.insert(n);
    }
    return output.size();
  }
}

