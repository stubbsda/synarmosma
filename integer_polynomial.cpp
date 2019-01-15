#include "integer_polynomial.h"

using namespace SYNARMOSMA;

extern Random RND;

template<class kind>
Integer_Polynomial<kind>::Integer_Polynomial()
{

}

template<class kind>
Integer_Polynomial<kind>::Integer_Polynomial(unsigned int n)
{
  degree = n;
  initialize();
}

template<class kind>
Integer_Polynomial<kind>::Integer_Polynomial(unsigned int n,unsigned int p)
{
  if (p > 0) {
    if (!NTL::ProbPrime(p)) throw std::invalid_argument("Non-zero field characteristic must be prime!");
  }
  characteristic = p;
  degree = n;
  initialize();
}
/*
namespace SYNARMOSMA {
  template<>
  void Integer_Polynomial<NTL::ZZ>::simplify()
  {
    unsigned int i,d = 0;
    for(i=0; i<=degree; ++i) {
      if (terms[i] == NTL::to_ZZ(0)) continue;
      d = i;
    }
    if (d < degree) {
      std::vector<NTL::ZZ> nterms;

      for(i=0; i<=d; ++i) {
        nterms.push_back(terms[i]);
      }
      terms = nterms;
      degree = d;
    }
    if (terms[degree] == NTL::to_ZZ(1)) monic = true;
    if (terms[0] == NTL::to_ZZ(0)) homogeneous = true;
  }
}

template<class kind>
void Integer_Polynomial<kind>::simplify()
{
  unsigned int i,d = 0;
  for(i=0; i<=degree; ++i) {
    if (terms[i] == kind(0)) continue;
    d = i;
  }
  if (d < degree) {
    std::vector<kind> nterms;

    for(i=0; i<=d; ++i) {
      nterms.push_back(terms[i]);
    }
    terms = nterms;
    degree = d;
  }
  if (terms[degree] == kind(1)) monic = true;
  if (terms[0] == kind(0)) homogeneous = true;
}
*/
template<class kind>
Integer_Polynomial<kind>::Integer_Polynomial(const std::vector<kind>& t)
{
  unsigned int i;

  degree = t.size() - 1;
  for(i=0; i<t.size(); ++i) {
    terms.push_back(t[i]);
  }

  simplify();
}

template<class kind>
Integer_Polynomial<kind>::Integer_Polynomial(const std::vector<kind>& t,unsigned int p)
{
  if (p > 0) {
    if (!NTL::ProbPrime(p)) throw std::invalid_argument("Non-zero field characteristic must be prime!");
  }
  unsigned int i;

  characteristic = p;
  degree = t.size() - 1;
  for(i=0; i<t.size(); ++i) {
    terms.push_back(t[i]);
  }
  irreducible = false;
  homogeneous = false;
  monic = false;
  simplify();
}

template<class kind>
Integer_Polynomial<kind>::Integer_Polynomial(const Integer_Polynomial& source)
{
  degree = source.degree;
  characteristic = source.characteristic;
  monic = source.monic;
  homogeneous = source.homogeneous;
  irreducible = source.irreducible;
  terms = source.terms;
}

template<class kind>
Integer_Polynomial<kind>::~Integer_Polynomial()
{

}

template<class kind>
void Integer_Polynomial<kind>::clear()
{
  terms.clear();
  characteristic = 0;
  degree = 0;
  homogeneous = false;
  irreducible = false;
  monic = false;
}

template<class kind>
void Integer_Polynomial<kind>::generate(unsigned int d)
{
  unsigned int chi = characteristic;
  clear();
  degree = d;
  characteristic = chi;
  initialize();
}

namespace SYNARMOSMA {
  template<>
  int Integer_Polynomial<NTL::ZZ>::write_terms(std::ofstream& s) const
  {
    unsigned int i;
    int count = 0;

    for(i=0; i<=degree; ++i) {
      count += write_ZZ(s,terms[i]);
    }
    return count;
  }
}

template<class kind>
int Integer_Polynomial<kind>::write_terms(std::ofstream& s) const
{
  unsigned int i;
  int count = 0;
  kind t;

  for(i=0; i<=degree; ++i) {
    t = terms[i];
    s.write((char*)(&t),sizeof(kind)); count += sizeof(kind);
  }
  return count;
}

template<class kind>
int Integer_Polynomial<kind>::serialize(std::ofstream& s) const
{
  int count = 0;

  s.write((char*)(&degree),sizeof(int)); count += sizeof(int);
  s.write((char*)(&characteristic),sizeof(int)); count += sizeof(int);
  s.write((char*)(&homogeneous),sizeof(bool)); count += sizeof(bool);
  s.write((char*)(&irreducible),sizeof(bool)); count += sizeof(bool);
  s.write((char*)(&monic),sizeof(bool)); count += sizeof(bool);
  count += write_terms(s);

  return count;
}

namespace SYNARMOSMA {
  template<>
  int Integer_Polynomial<NTL::ZZ>::read_terms(std::ifstream& s)
  {
    unsigned int i;
    int count = 0;
    NTL::ZZ t;

    for(i=0; i<=degree; ++i) {
      count += read_ZZ(s,t);
      terms.push_back(t);
    }
    return count;
  }
}

template<class kind>
int Integer_Polynomial<kind>::read_terms(std::ifstream& s)
{
  unsigned int i;
  int count = 0;
  kind t;

  for(i=0; i<=degree; ++i) {
    s.read((char*)(&t),sizeof(kind)); count += sizeof(kind);
    terms.push_back(t);
  }

  return count;
}

template<class kind>
int Integer_Polynomial<kind>::deserialize(std::ifstream& s)
{
  int count = 0;

  clear();

  s.read((char*)(&degree),sizeof(int)); count += sizeof(int);
  s.read((char*)(&characteristic),sizeof(int)); count += sizeof(int);
  s.read((char*)(&homogeneous),sizeof(bool)); count += sizeof(bool);
  s.read((char*)(&irreducible),sizeof(bool)); count += sizeof(bool);
  s.read((char*)(&monic),sizeof(bool)); count += sizeof(bool);
  count += read_terms(s);

  return count;
}

template<class kind>
void Integer_Polynomial<kind>::initialize()
{
  unsigned int i;

  for(i=0; i<=degree; ++i) {
    terms.push_back(Integer_Polynomial<kind>::unity);
  }

  irreducible = true;
  simplify();
}

template<class kind>
Integer_Polynomial<kind>& Integer_Polynomial<kind>::operator =(const Integer_Polynomial<kind>& source)
{
  if (this == &source) return *this;

  degree = source.degree;
  characteristic = source.characteristic;
  monic = source.monic;
  homogeneous = source.homogeneous;
  irreducible = source.irreducible;
  terms = source.terms;

  return *this;
}

template<class kind>
Integer_Polynomial<kind>& Integer_Polynomial<kind>::operator -(const Integer_Polynomial<kind>& source)
{
  unsigned int i;

  degree = source.degree;
  characteristic = source.characteristic;
  monic = source.monic;
  homogeneous = source.homogeneous;
  irreducible = source.irreducible;
  terms = source.terms;  
  for(i=0; i<=degree; ++i) {
    terms[i] = neg1*terms[i];
  }

  return *this;
}

template<class kind>
kind Integer_Polynomial<kind>::get_value(unsigned int n) const
{
  assert(n <= degree);
  return terms[n];
}

template<class kind>
void Integer_Polynomial<kind>::set_value(kind x,unsigned int n)
{
  if (n <= degree) {
    terms[n] = x;
  }
  else {
    unsigned int i;
    for(i=1+degree; i<n; ++i) {
      terms.push_back(Integer_Polynomial<kind>::zero);
    }
    terms.push_back(x);
    degree = n;
  }
  simplify();
}

namespace SYNARMOSMA {
  template<>
  NTL::ZZ Integer_Polynomial<NTL::ZZ>::evaluate(NTL::ZZ x)
  {
    if (characteristic > 0) {
      if (x >= characteristic) throw std::invalid_argument("Argument exceeds field characteristic!");  
    }
    int i;
    NTL::ZZ y = NTL::to_ZZ(0);
    for(i=degree; i>=0; --i) {
      y = y*x + terms[i];
    }
    if (characteristic > 1) y = y % characteristic;    
    return y;
  }
}

template<class kind>
kind Integer_Polynomial<kind>::evaluate(kind x)
{
  if (characteristic > 0) {
    if (x >= characteristic) throw std::invalid_argument("Argument exceeds field characteristic!");
  }
  int i;
  kind y = Integer_Polynomial<kind>::zero;
  // Use Horner's method to speed evaluation of the polynomial
  for(i=degree; i>0; --i) {
    y = y*x + terms[i];
  }
  if (characteristic > 1) y = y % characteristic;
  return y;
}

namespace SYNARMOSMA {
  template<>
  Integer_Polynomial<NTL::ZZ> Integer_Polynomial<NTL::ZZ>::derivative() const
  {
    unsigned int i;
    Integer_Polynomial<NTL::ZZ> output(degree-1);
    for(i=0; i<degree-1; ++i) {
      output.set_value(NTL::to_ZZ(1+i)*terms[i+1],i);
    }
    output.simplify();
    return output;
  }
}

template<class kind>
Integer_Polynomial<kind> Integer_Polynomial<kind>::derivative() const
{
  unsigned int i;
  Integer_Polynomial<kind> output(degree-1);
  for(i=0; i<degree-1; ++i) {
    output.set_value(kind(1+i)*terms[i+1],i);
  }
  output.simplify();
  return output;
}

namespace SYNARMOSMA {
  template<>
  Integer_Polynomial<int> Integer_Polynomial<NTL::ZZ>::reduce(unsigned int p)
  {
    if (!NTL::ProbPrime(p)) throw std::invalid_argument("Integer_Polynomial must be reduced over a prime characteristic!");

    unsigned int i;
    long temp;
    NTL::ZZ q,base = NTL::to_ZZ(p);
    std::vector<int> nterms;

    for(i=0; i<=degree; ++i) {
      q = terms[i] % base;
      NTL::conv(temp,q);
      nterms.push_back(temp);
    }
    Integer_Polynomial<int> output(nterms,p);
    return output;
  }
}

template<class kind>
Integer_Polynomial<int> Integer_Polynomial<kind>::reduce(unsigned int p)
{
  if (!NTL::ProbPrime(p)) throw std::invalid_argument("Integer_Polynomial must be reduced over a prime characteristic!");
  
  unsigned int i;
  std::vector<int> nterms;

  for(i=0; i<=degree; ++i) {
    nterms.push_back(convert(terms[i],p));
  }
  Integer_Polynomial<int> output(nterms,p);
  return output;
}



