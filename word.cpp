#include "word.h"

using namespace SYNARMOSMA;

extern Random RND;

Word::Word()
{

}

Word::Word(int p)
{
  assert(p > 0);
  NL = (unsigned) p;
  unsigned int n = 4 + RND.irandom(20);
  initialize(n);
}

Word::Word(int p,int n)
{
  assert(p > 0 && n > 0);
  NL = (unsigned) p;
  initialize((unsigned) n);
}

Word::Word(int p,const std::string& w)
{
  assert(p > 0);
  unsigned int i,j,n;
  int e;
  char alphabet[] = {'a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z'};

  NL = (unsigned) p;
  for(i=0; i<w.length(); ++i) {
    assert(std::isalpha(w[i]));
    for(j=0; j<26; ++j) {
      if (std::tolower(w[i]) == alphabet[j]) {
        n = j;
        break;
      }
    }
    e = (std::isupper(w[i])) ? -1 : 1;
    content.push_back(std::pair<unsigned int,int>(n,e));
  }
}

Word::Word(int p,int n,int m)
{
  assert(p > 0 && n > 0);

  NL = (unsigned) p;
#ifdef DEBUG
  assert(n < p);
#endif
  std::pair<unsigned int,int> doublet((unsigned) n,m);
  content.push_back(doublet);
}

Word::Word(const Word& source)
{
  content = source.content;
  NL = source.NL;
}

Word& Word::operator =(const Word& source)
{
  if (this == &source) return *this;

  content = source.content;
  NL = source.NL;

  return *this;
}

Word::~Word()
{

}

bool Word::legal() const
{
  std::vector<std::pair<unsigned int,int> >::const_iterator it;

  for(it=content.begin(); it!=content.end(); ++it) {
    if (it->first >= NL) return false;
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
    if ((std::abs(content[0].second) == 1) && (std::abs(content[1].second) == 1)) return true;
  }
  return false;
}

bool Word::homogeneous() const
{
  // This is needed for the algorithm to test words in a braid group for their 
  // equivalence
  unsigned int prime = 10000;
  bool ppow = false,npow = false;
  std::vector<std::pair<unsigned int,int> >::const_iterator it;

  for(it=content.begin(); it!=content.end(); ++it) {
    if (it->first < prime) prime = it->first;
  }
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

void Word::permute(int n,Word& w) const
{
  // This method will create the word formed by the cyclic
  // permutation of *this, that is
  // a_1*a_2*...*a_k -> a_n*a_{n+1}*...*a_k*a_1*a_2*...*a_{n-1}
  assert(n >= 0);

  w.clear();
#ifdef DEBUG
  assert(n < (signed) content.size());
#endif
  if (n == 0) {
    w = Word(*this);
    return;
  }
  for(int i=n; i<(signed) content.size(); ++i) {
    w.content.push_back(content[i]);
  }
  for(int i=0; i<n; ++i) {
    w.content.push_back(content[i]);
  }
}

void Word::initialize(unsigned int n)
{
  unsigned int i;
  std::pair<unsigned int,int> doublet;

  content.clear();
  for(i=0; i<n; ++i) {
    doublet.first = (unsigned) RND.irandom(NL);
    doublet.second = RND.irandom(1,10);
    if (RND.drandom() < 0.5) doublet.second *= -1;
    content.push_back(doublet);
  }
}

void Word::initialize(unsigned int p,unsigned int n,int m)
{
  NL = p;
#ifdef DEBUG
  assert(n < p);
#endif
  std::pair<unsigned int,int> doublet(n,m);
  content.push_back(doublet);
}

void Word::initialize(const std::vector<unsigned int>& base,const std::vector<int>& exponent)
{
  unsigned int i;
  std::pair<unsigned int,int> doublet;
#ifdef DEBUG
  assert(base.size() == exponent.size());
#endif
  content.clear();
  for(i=0; i<base.size(); ++i) {
#ifdef DEBUG
    assert(base[i] < NL);
#endif
    doublet.first = base[i];
    doublet.second = exponent[i];
    content.push_back(doublet);
  }
}

Word Word::operator !() const
{
  return (this->invert());
}

Word Word::invert() const
{  
  int i,n = (signed) content.size();
  Word output(NL,0);
  std::pair<unsigned int,int> doublet;

  for(i=n-1; i>=0; --i) {
    doublet.first = content[i].first;
    doublet.second = -content[i].second;
    output.content.push_back(doublet);
  }
  return output;
}

Word Word::swap(int p,int q,bool inv) const
{
  assert(p >= 0 && q >= 0);
  unsigned int i,n = content.size();
  Word output(0);
  std::pair<unsigned int,int> doublet;

  for(i=0; i<n; ++i) {
    doublet.first = content[i].first;
    doublet.second = content[i].second;
    if (doublet.first == (unsigned) q) {
      doublet.first = (unsigned) p;
      if (inv) doublet.second *= -1;
    }
    else if (doublet.first > (unsigned) q) {
      doublet.first -= 1;
    }
    output.content.push_back(doublet);
  }
  return output;
}

Word Word::mutate() const
{
  Word output(0);
  std::pair<unsigned int,int> doublet;
  int n = RND.irandom(content.size());

  output.content = content;

  doublet.first = RND.irandom(NL);
  doublet.second = RND.irandom(1,10);
  if (RND.drandom() < 0.5) doublet.second = -doublet.second;
  output.content[n] = doublet;

  return output;
}

int Word::serialize(std::ofstream& s) const
{
  unsigned int i,j,n = content.size();
  int k,count = 0;

  s.write((char*)(&NL),sizeof(int)); count += sizeof(int);
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

  s.read((char*)(&NL),sizeof(int)); count += sizeof(int);
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

void Word::free_reduce() 
{
  int i,rpoint = -1,L = (signed) content.size();
  unsigned int n,m;
  std::vector<std::pair<unsigned int,int> > vx;

  do {
    rpoint = -1;
    L = (signed) content.size();
    for(i=0; i<L-1; ++i) {
      n = content[i].first;
      m = content[i+1].first;
      if (n == m) {
        if (content[i].second == -content[i+1].second) {
          rpoint = i;
          break;
        }
      }
    }
    if (rpoint == -1) break;
    vx.clear();
    for(i=0; i<L; ++i) {
      if (i == rpoint || i == (rpoint+1)) continue;
      vx.push_back(content[i]);
    }
    content = vx;
  } while(true);
}

Word Word::reduce(int M,const std::set<unsigned int>& trivial_generators,const unsigned int* offset) const
{
  Word output(M,0);
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
  const unsigned int n = content.size();
  std::pair<unsigned int,int> doublet;
  Word output(NL,0);

  if (content.size() < 2) {
    for(i=0; i<content.size(); ++i) {
      output.content.push_back(content[i]);
    }
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
  return output;
}

namespace SYNARMOSMA {
  bool operator ==(const Word& w1,const Word& w2)
  {
    unsigned int i;

    if (w1.content.size() == w2.content.size()) {
      bool good = true;
      for(i=0; i<w1.content.size(); ++i) {
        if (w1.content[i] != w2.content[i]) {
          good = false;
          break;
        }
      }
      if (good) return true;
    }
    return false;
  }

  bool operator !=(const Word& w1,const Word& w2)
  {
    if (w1 == w2) return false;
    return true;
  }

  Word operator *(const Word& w1,const Word& w2)
  {
    unsigned int i;
    unsigned int n = (w1.NL > w2.NL) ? w1.NL : w2.NL;
    Word output(n,0);
 
    for(i=0; i<w1.content.size(); ++i) {
      output.content.push_back(w1.content[i]);
    }
    for(i=0; i<w2.content.size(); ++i) {
      output.content.push_back(w2.content[i]);
    }
    return output;
  }

  std::ostream& operator <<(std::ostream& os,const Word& source)
  {
    unsigned int i;
    std::pair<unsigned int,int> doublet;

    if (source.content.empty()) return os;

    for(i=0; i<source.content.size()-1; ++i) {
      doublet = source.content[i];
      os << "x[" << 1 + doublet.first << "]^(" << doublet.second << ")*";
    }
    doublet = source.content[source.content.size()-1];
    os << "x[" << 1 + doublet.first << "]^(" << doublet.second << ")";
    return os;
  }
}

