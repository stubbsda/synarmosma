#include "functional_equation.h"

using namespace SYNARMOSMA;

template<class kind>
Functional_Equation<kind>::Functional_Equation()
{

}

template<class kind>
Functional_Equation<kind>::Functional_Equation(unsigned int n)
{
  if (n == 0) throw std::invalid_argument("The degree of the functional equation must be greater than zero!");
  initialize(n);
}

template<class kind>
Functional_Equation<kind>::Functional_Equation(const std::string& filename)
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
      if (bk.size() != 2) throw std::runtime_error("Error parsing the functional equation input file!");
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
  void Functional_Equation<NTL::ZZ>::analyze_file(std::vector<std::string>& alpha,std::vector<std::string>& beta,std::vector<std::string>& exponents)
  {
    unsigned int i,j,degree;
    char c;
    std::string store;
    Integer_Polynomial<NTL::ZZ> p1,p2;
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
      p1 = Integer_Polynomial<NTL::ZZ>(vx);

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
      p2 = Integer_Polynomial<NTL::ZZ>(vx);

      vx.clear();
      store.clear();
      std::tuple<Integer_Polynomial<NTL::ZZ>,Integer_Polynomial<NTL::ZZ>,unsigned int> trio(p1,p2,degree);
      terms.push_back(trio);
    }
    simplify();
    if (!consistent()) throw std::invalid_argument("The functional equation input is not consistent with this class!");
  }
}

template<class kind>
void Functional_Equation<kind>::analyze_file(std::vector<std::string>& alpha,std::vector<std::string>& beta,std::vector<std::string>& exponents)
{
  unsigned int i,j,degree;
  char c;
  std::string store;
  Integer_Polynomial<kind> p1,p2;
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
    p1 = Integer_Polynomial<kind>(vx);

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
    p2 = Integer_Polynomial<kind>(vx);

    vx.clear();
    store.clear();
    std::tuple<Integer_Polynomial<kind>,Integer_Polynomial<kind>,unsigned int> trio(p1,p2,degree);
    terms.push_back(trio);
  }
  simplify();
  if (!consistent()) throw std::invalid_argument("The functional equation input is not consistent with this class!");
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
  Integer_Polynomial<kind> p;
  
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
  Integer_Polynomial<kind> p1,p2;

  clear();

  s.read((char*)(&linear),sizeof(bool)); count += sizeof(bool);
  s.read((char*)(&homogeneous),sizeof(bool)); count += sizeof(bool);
  s.read((char*)(&n),sizeof(int)); count += sizeof(int);
  for(i=0; i<n; ++i) {
    s.read((char*)(&j),sizeof(int)); count += sizeof(int);
    count += p1.deserialize(s);
    count += p2.deserialize(s);
    terms.push_back(std::tuple<Integer_Polynomial<kind>,Integer_Polynomial<kind>,unsigned int>(p1,p2,j));
  }
  count += remainder.deserialize(s);

  return count;
}

template<class kind>
void Functional_Equation<kind>::initialize(unsigned int n)
{
  unsigned int i,d = (5 < (1 + 2*n)) ? 5 : 1 + 2*n;
  Integer_Polynomial<kind> q,p;

  for(i=1; i<n; ++i) {
    q.generate(d - 1);
    p.generate(d);
    terms.push_back(std::tuple<Integer_Polynomial<kind>,Integer_Polynomial<kind>,unsigned int>(q,p,i));
  }
  do {
    q.generate(d - 1);
  } while(q.is_null());
  p.generate(d);
  terms.push_back(std::tuple<Integer_Polynomial<kind>,Integer_Polynomial<kind>,unsigned int>(q,p,n));

  q.generate(d - 1);
  remainder = q;

  simplify();
  if (!consistent()) throw std::runtime_error("The randomly generated functional equation is inconsistent!");
}

template<class kind>
bool Functional_Equation<kind>::consistent() const
{
  unsigned int i,d,n = terms.size();
  std::set<unsigned int> dset;

  // Begin by checking that all of the terms are distinct in terms of their degree:
  for(i=0; i<n; ++i) {
    d = std::get<2>(terms[i]);
    if (dset.count(d) > 0) return false;
    dset.insert(d);
  }
  return true;
}

template<class kind>
bool Functional_Equation<kind>::simplify()
{
  unsigned int i,d = 0,n = terms.size();
  bool output = false;
  Integer_Polynomial<kind> q;
  std::vector<std::tuple<Integer_Polynomial<kind>,Integer_Polynomial<kind>,unsigned int> > nterms;

  if (!homogeneous) {
    if (remainder.is_null()) {
      homogeneous = true;
      output = true;
    }
  }
  for(i=0; i<n; ++i) {
    q = std::get<0>(terms[i]);
    if (q.is_null()) continue;
    if (std::get<2>(terms[i]) > d) d = std::get<2>(terms[i]);
    nterms.push_back(terms[i]);
  }
  if (d == 1) {
    if (!linear) {
      linear = true;
      output = true;
    }
  }
  if (nterms.size() == terms.size()) return output;
  terms = nterms;
  return true;  
}

namespace SYNARMOSMA
{
  template<>
  Variety<unsigned int> Functional_Equation<NTL::ZZ>::reduce(unsigned int p)
  {
    if (!NTL::ProbPrime(p)) throw std::invalid_argument("Functional equation must be reduced over a prime!");
    unsigned int i,j,in1;
    long q;
    NTL::ZZ z;
    std::pair<unsigned int,unsigned int> duo;
    std::tuple<Integer_Polynomial<NTL::ZZ>,Integer_Polynomial<NTL::ZZ>,unsigned int> trio;
    Variety<unsigned int> output(p,p);
    Integer_Polynomial<NTL::ZZ> py;
    Monomial<unsigned int> term;

    output.clear();

    for(i=0; i<p; ++i) {
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
Variety<unsigned int> Functional_Equation<kind>::reduce(unsigned int p)
{
  if (!NTL::ProbPrime(p)) throw std::invalid_argument("Functional equation must be reduced over a prime!");

  unsigned int i,j,in1;
  std::pair<unsigned int,unsigned int> duo;
  std::tuple<Integer_Polynomial<kind>,Integer_Polynomial<kind>,unsigned int> trio;
  Variety<unsigned int> output(p,p);
  Integer_Polynomial<kind> py;
  Monomial<unsigned int> term;

  output.clear();

  for(i=0; i<p; ++i) {
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
