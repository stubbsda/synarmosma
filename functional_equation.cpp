#include "functional_equation.h"

using namespace SYNARMOSMA;

extern Random RND;

template<class kind>
Functional_Equation<kind>::Functional_Equation()
{
  unsigned int n = 5 + RND.irandom(6);

  linear = false;
  homogeneous = false;
  float alpha = float(RND.drandom());
  if (alpha < 0.5) linear = true;
  alpha = float(RND.drandom());
  if (alpha < 0.5) homogeneous = true;
  initialize(n);
}

template<class kind>
Functional_Equation<kind>::Functional_Equation(int n)
{
  assert(n > 0);
  initialize(n);
}

template<class kind>
Functional_Equation<kind>::Functional_Equation(const char* filename)
{
  unsigned int i;
  bool first = true;
  std::string line,store;
  std::vector<std::string> alpha,beta,exponents;
  std::vector<unsigned int> bk;

  std::ifstream s;
  s.exceptions(std::ifstream::badbit);
  try {
    s.open(filename,std::ios::in);

    // Loop through all lines in the parameter file
    while(std::getline(s,line)) {
      // If it's an empty line, continue
      if (line.empty()) continue;
      // If the line begins with a #, ignore it
      if (line[0] == '#') continue;
      if (first) {
        // We need to determine the field of this equation 
        trim(line);
        first = false;
        continue;
      }
      // Find the position of the colons
      for(i=0; i<line.length(); ++i) {
        if (line[i] == ':') bk.push_back(i);
      }
#ifdef DEBUG
      assert(bk.size() == 2);
#endif
      store = line.substr(0,bk[0]);
      trim(store);
      alpha.push_back(store);
      store = line.substr(bk[0]+1,bk[1]-1);
      trim(store);
      beta.push_back(store);
      store = line.substr(bk[1]+1,line.length());
      trim(store);
      exponents.push_back(store);
      bk.clear();
    }
  }
  catch (const std::ifstream::failure& e) {
    std::cout << "Error in opening or reading the " << filename << " file!" << std::endl;
  }
  s.close();

  analyze_file(alpha,beta,exponents);
}

namespace SYNARMOSMA {
  template<>
  void Functional_Equation<Rational>::analyze_file(std::vector<std::string>& alpha,std::vector<std::string>& beta,std::vector<std::string>& exponents)
  {
    Polynomial<Rational> p1,p2;
    std::vector<Rational> vx;
    Rational coefficient;
    unsigned int i,j,degree;
    char c;
    std::string s1,s2;
    bool numerator;

    for(i=0; i<alpha.size(); ++i) {
      degree = boost::lexical_cast<unsigned int>(exponents[i]);
      numerator = true;
      for(j=0; j<alpha[i].length(); ++j) {
        c = alpha[i][j];
        if (c == '(') continue;
        if (c == ',' || c == ')') {
          coefficient = Rational(boost::lexical_cast<signed int>(s1),boost::lexical_cast<signed int>(s2));
          vx.push_back(coefficient);
          s1.clear();
          s2.clear();
          continue;
        }
        if (c == '/') {
          numerator = false;
          continue;
        }
        if (numerator) {
          s1.push_back(c);
        }
        else {
          s2.push_back(c);
        }
      }
      p1 = Polynomial<Rational>(vx);

      vx.clear();
      s1.clear();
      s2.clear();
      if (degree == 0) {
        remainder = p1;
        continue;
      }
      numerator = true;
      for(j=0; j<beta[i].length(); ++j) {
        c = beta[i][j];
        if (c == '(') continue;
        if (c == ',' || c == ')') {
          coefficient = Rational(boost::lexical_cast<signed int>(s1),boost::lexical_cast<signed int>(s2));
          vx.push_back(coefficient);
          s1.clear();
          s2.clear();
          continue;
        }
        if (c == '/') {
          numerator = false;
          continue;
        }
        if (numerator) {
          s1.push_back(c);
        }
        else {
          s2.push_back(c);
        }
      }
      p2 = Polynomial<Rational>(vx);

      std::tuple<Polynomial<Rational>,Polynomial<Rational>,unsigned int> trio(p1,p2,degree);
      terms.push_back(trio);
      vx.clear();
      s1.clear();
      s2.clear();
    }
  }

  template<>
  void Functional_Equation<NTL::ZZ>::analyze_file(std::vector<std::string>& alpha,std::vector<std::string>& beta,std::vector<std::string>& exponents)
  {
    unsigned int i,j,degree;
    char c;
    std::string store;
    Polynomial<NTL::ZZ> p1,p2;
    std::vector<NTL::ZZ> vx;

    for(i=0; i<alpha.size(); ++i) {
      degree = boost::lexical_cast<unsigned int>(exponents[i]);
      for(j=0; j<alpha[i].length(); ++j) {
        c = alpha[i][j];
        if (c == '(') continue;
        if (c == ',' || c == ')') {
          vx.push_back(NTL::to_ZZ(boost::lexical_cast<signed int>(store)));
          store.clear();
          continue;
        }
        store.push_back(c);
      }
      p1 = Polynomial<NTL::ZZ>(vx);

      vx.clear();
      store.clear();
      if (degree == 0) {
        remainder = p1;
        continue;
      }
      for(j=0; j<beta[i].length(); ++j) {
        c = beta[i][j];
        if (c == '(') continue;
        if (c == ',' || c == ')') {
          vx.push_back(NTL::to_ZZ(boost::lexical_cast<signed int>(store)));
          store.clear();
          continue;
        }
        store.push_back(c);
      }
      p2 = Polynomial<NTL::ZZ>(vx);

      vx.clear();
      store.clear();
      std::tuple<Polynomial<NTL::ZZ>,Polynomial<NTL::ZZ>,unsigned int> trio(p1,p2,degree);
      terms.push_back(trio);
    }
  }
}

template<class kind>
void Functional_Equation<kind>::analyze_file(std::vector<std::string>& alpha,std::vector<std::string>& beta,std::vector<std::string>& exponents)
{
  unsigned int i,j,degree;
  char c;
  std::string store;
  Polynomial<kind> p1,p2;
  std::vector<kind> vx;

  for(i=0; i<alpha.size(); ++i) {
    degree = boost::lexical_cast<unsigned int>(exponents[i]);
    for(j=0; j<alpha[i].length(); ++j) {
      c = alpha[i][j];
      if (c == '(') continue;
      if (c == ',' || c == ')') {
        vx.push_back(boost::lexical_cast<signed int>(store));
        store.clear();
        continue;
      }
      store.push_back(c);
    }
    p1 = Polynomial<kind>(vx);

    vx.clear();
    store.clear();
    if (degree == 0) {
      remainder = p1;
      continue;
    }
    for(j=0; j<beta[i].length(); ++j) {
      c = beta[i][j];
      if (c == '(') continue;
      if (c == ',' || c == ')') {
        vx.push_back(boost::lexical_cast<signed int>(store));
        store.clear();
        continue;
      }
      store.push_back(c);
    }
    p2 = Polynomial<kind>(vx);

    vx.clear();
    store.clear();
    std::tuple<Polynomial<kind>,Polynomial<kind>,unsigned int> trio(p1,p2,degree);
    terms.push_back(trio);
  }
}

template<class kind>
Functional_Equation<kind>::Functional_Equation(const Functional_Equation& source)
{
  linear = source.linear;
  homogeneous = source.homogeneous;
  terms = source.terms;
  remainder = source.remainder;
}

template<class kind>
Functional_Equation<kind>& Functional_Equation<kind>::operator =(const Functional_Equation<kind>& source)
{
  if (this == &source) return *this;
  linear = source.linear;
  homogeneous = source.homogeneous;
  terms = source.terms;
  remainder = source.remainder;
  return *this;
}

template<class kind>
Functional_Equation<kind>::~Functional_Equation()
{

}

template<class kind>
void Functional_Equation<kind>::clear()
{
  terms.clear();
  remainder.clear();
  linear = false;
  homogeneous = false;
}

template<class kind>
int Functional_Equation<kind>::serialize(std::ofstream& s) const
{
  unsigned int i,j,n;
  int count = 0;
  Polynomial<kind> p;
  
  s.write((char*)(&linear),sizeof(bool)); count += sizeof(bool);
  s.write((char*)(&homogeneous),sizeof(bool)); count += sizeof(bool);
  n = terms.size();
  s.write((char*)(&n),sizeof(int)); count += sizeof(int);
  for(i=0; i<n; ++i) {
    j = std::get<2>(terms[i]);
    s.write((char*)(&j),sizeof(int)); count += sizeof(int);
    p = std::get<0>(terms[i]);
    count += p.serialize(s);
    p = std::get<1>(terms[i]);
    count += p.serialize(s);
  }
  count += remainder.serialize(s);
  
  return count;
}

template<class kind>
int Functional_Equation<kind>::deserialize(std::ifstream& s)
{
  unsigned int i,j,n;
  int count = 0;
  Polynomial<kind> p1,p2;

  clear();

  s.read((char*)(&linear),sizeof(bool)); count += sizeof(bool);
  s.read((char*)(&homogeneous),sizeof(bool)); count += sizeof(bool);
  s.read((char*)(&n),sizeof(int)); count += sizeof(int);
  for(i=0; i<n; ++i) {
    s.read((char*)(&j),sizeof(int)); count += sizeof(int);
    count += p1.deserialize(s);
    count += p2.deserialize(s);
    terms.push_back(std::tuple<Polynomial<kind>,Polynomial<kind>,unsigned int>(p1,p2,j));
  }
  count += remainder.deserialize(s);

  return count;
}

template<class kind>
void Functional_Equation<kind>::initialize(int n)
{
  assert(n > 0);
  int i;
  unsigned int m = 1;
  for(i=1; i<n; ++i) {
    Polynomial<kind> alpha,beta;
    if (!linear) m = 1 + RND.irandom(8);
    std::tuple<Polynomial<kind>,Polynomial<kind>,unsigned int> triple(alpha,beta,m);
    terms.push_back(triple);
  }
  if (!homogeneous) {
    Polynomial<kind> alpha;
    remainder = alpha;
  }
}

namespace SYNARMOSMA
{
  template<>
  Variety<unsigned int> Functional_Equation<NTL::ZZ>::reduce(int p)
  {
    assert(p > 0);
    if (!NTL::ProbPrime(p)) throw std::invalid_argument("Functional equation must be reduced over a prime!");
    unsigned int j,in1;
    long q;
    NTL::ZZ z;
    std::pair<unsigned int,unsigned int> duo;
    std::tuple<Polynomial<NTL::ZZ>,Polynomial<NTL::ZZ>,unsigned int> trio;
    Variety<unsigned int> output(p,p);
    Polynomial<NTL::ZZ> py;
    Monomial<unsigned int> term;

    output.clear();

    for(int i=0; i<p; ++i) {
      z = NTL::to_ZZ(long(i));
      for(j=0; j<terms.size(); ++j) {
        trio = terms[j];
        // First convert alpha(p) to an element of GF(p)
        py = std::get<0>(trio);
        NTL::conv(q,py.evaluate(z));
        in1 = (unsigned int) q;
        in1 = in1 % p;
        if (in1 == 0) continue;
        term.coefficient = in1;
        // Next we need to convert beta(p) to an element of GF(p)
        py = std::get<1>(trio);
        NTL::conv(q,py.evaluate(z));
        in1 = (unsigned int) q;
        in1 = in1 % p;
        duo.first = in1;
        // Finally we have to find out the exponent 
        duo.second = std::get<2>(trio);
        term.exponents.push_back(duo);
        output.add_term(i,term);
        term.exponents.clear();
      }
      NTL::conv(q,remainder.evaluate(z));
      in1 = (unsigned int) q;
      in1 = in1 % p;
      output.set_remainder_value(i,in1);
    }
    output.elaborate();
    return output;
  }
}

template<class kind>
Variety<unsigned int> Functional_Equation<kind>::reduce(int p)
{
  assert(p > 0);
  unsigned int j,in1;
  std::pair<unsigned int,unsigned int> duo;
  std::tuple<Polynomial<kind>,Polynomial<kind>,unsigned int> trio;
  Variety<unsigned int> output(p,p);
  Polynomial<kind> py;
  Monomial<unsigned int> term;

  output.clear();

  for(int i=0; i<p; ++i) {
    for(j=0; j<terms.size(); ++j) {
      trio = terms[j];
      // First convert alpha(p) to an element of GF(p)
      py = std::get<0>(trio);
      in1 = convert(py.evaluate(i),p);
      if (in1 == 0) continue;
      term.coefficient = in1;
      // Next we need to convert beta(p) to an element of GF(p)
      py = std::get<1>(trio);
      duo.first = convert(py.evaluate(i),p);
      // Finally we have to find out the exponent 
      duo.second = std::get<2>(trio);
      term.exponents.push_back(duo);
      output.add_term(i,term);
      term.exponents.clear();
    }
    in1 = convert(remainder.evaluate(i),p);
    output.set_remainder_value(i,in1);
  }
  output.elaborate();
  return output;
}
