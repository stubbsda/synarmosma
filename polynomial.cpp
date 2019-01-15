#include "polynomial.h"

using namespace SYNARMOSMA;

extern Random RND;

template<class kind>
Polynomial<kind>::Polynomial()
{

}

template<class kind>
Polynomial<kind>::Polynomial(unsigned int n)
{
  degree = n;
  initialize();
}

template<class kind>
Polynomial<kind>::Polynomial(const std::vector<kind>& t)
{
  unsigned int i;

  degree = t.size() - 1;
  for(i=0; i<t.size(); ++i) {
    terms.push_back(t[i]);
  }

  simplify();
}

template<class kind>
Polynomial<kind>::Polynomial(const Polynomial& source)
{
  degree = source.degree;
  monic = source.monic;
  homogeneous = source.homogeneous;
  irreducible = source.irreducible;
  terms = source.terms;
}

template<class kind>
Polynomial<kind>::~Polynomial()
{

}

template<class kind>
void Polynomial<kind>::clear()
{
  terms.clear();
  degree = 0;
  homogeneous = false;
  irreducible = false;
  monic = false;
}

template<class kind>
void Polynomial<kind>::generate(unsigned int d)
{
  clear();

  degree = d;
  initialize();
}

namespace SYNARMOSMA {
  template<>
  int Polynomial<Rational>::write_terms(std::ofstream& s) const
  {
    unsigned int i;
    int count = 0;

    for(i=0; i<=degree; ++i) {
      count += write_ZZ(s,terms[i].get_numerator());
      count += write_ZZ(s,terms[i].get_denominator());
    }
    return count;
  }
}

template<class kind>
int Polynomial<kind>::write_terms(std::ofstream& s) const
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
int Polynomial<kind>::serialize(std::ofstream& s) const
{
  int count = 0;

  s.write((char*)(&degree),sizeof(int)); count += sizeof(int);
  s.write((char*)(&homogeneous),sizeof(bool)); count += sizeof(bool);
  s.write((char*)(&irreducible),sizeof(bool)); count += sizeof(bool);
  s.write((char*)(&monic),sizeof(bool)); count += sizeof(bool);
  count += write_terms(s);

  return count;
}

namespace SYNARMOSMA {
  template<>
  int Polynomial<Rational>::read_terms(std::ifstream& s)
  {
    unsigned int i;
    int count = 0;
    NTL::ZZ t1,t2;

    for(i=0; i<=degree; ++i) {
      count += read_ZZ(s,t1);
      count += read_ZZ(s,t2);
      terms.push_back(Rational(t1,t2));
    }
    return count;
  }
}

template<class kind>
int Polynomial<kind>::read_terms(std::ifstream& s)
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
int Polynomial<kind>::deserialize(std::ifstream& s)
{
  int count = 0;

  clear();

  s.read((char*)(&degree),sizeof(int)); count += sizeof(int);
  s.read((char*)(&homogeneous),sizeof(bool)); count += sizeof(bool);
  s.read((char*)(&irreducible),sizeof(bool)); count += sizeof(bool);
  s.read((char*)(&monic),sizeof(bool)); count += sizeof(bool);
  count += read_terms(s);

  return count;
}

template<class kind>
void Polynomial<kind>::initialize()
{
  unsigned int i;
  kind test;
  bool flag = true;
  const kind zero = kind(0);

  for(i=0; i<degree; ++i) {
    terms.push_back(kind(RND.irandom(-25,25)));
  }
  do {
    test = kind(RND.irandom(-25,25));
    if (test != zero) flag = false;
  } while(flag);
  terms.push_back(test);

  irreducible = true;
  simplify();
}

template<class kind>
Polynomial<kind>& Polynomial<kind>::operator =(const Polynomial<kind>& source)
{
  if (this == &source) return *this;

  degree = source.degree;
  monic = source.monic;
  homogeneous = source.homogeneous;
  irreducible = source.irreducible;
  terms = source.terms;

  return *this;
}

template<class kind>
Polynomial<kind>& Polynomial<kind>::operator -(const Polynomial<kind>& source)
{
  unsigned int i;

  degree = source.degree;
  monic = source.monic;
  homogeneous = source.homogeneous;
  irreducible = source.irreducible;
  terms = source.terms;  
  for(i=0; i<=degree; ++i) {
    terms[i] = kind(-1)*terms[i];
  }

  return *this;
}

template<class kind>
kind Polynomial<kind>::get_value(unsigned int n) const
{
  assert(n <= degree);
  return terms[n];
}

template<class kind>
void Polynomial<kind>::set_value(kind x,unsigned int n)
{
  if (n <= degree) {
    terms[n] = x;
  }
  else {
    unsigned int i;
    for(i=1+degree; i<n; ++i) {
      terms.push_back(kind(0));
    }
    terms.push_back(x);
    degree = n;
  }
  simplify();
}

template<class kind>
kind Polynomial<kind>::evaluate(kind x)
{
  int i;
  kind y = kind(0);
  // Use Horner's method to speed evaluation of the polynomial
  for(i=degree; i>0; --i) {
    y = y*x + terms[i];
  }
  return y;
}

template<class kind>
Polynomial<kind> Polynomial<kind>::derivative() const
{
  unsigned int i;
  Polynomial<kind> output(degree-1);
  for(i=0; i<degree-1; ++i) {
    output.set_value(kind(1+i)*terms[i+1],i);
  }
  output.simplify();
  return output;
}

template<class kind>
Integer_Polynomial<unsigned int> Polynomial<kind>::reduce(unsigned int p)
{
  if (!NTL::ProbPrime(p)) throw std::invalid_argument("Polynomial must be reduced over a prime characteristic!");
  
  unsigned int i;
  std::vector<unsigned int> nterms;

  for(i=0; i<=degree; ++i) {
    nterms.push_back(convert(terms[i],p));
  }
  Integer_Polynomial<unsigned int> output(nterms,p);
  return output;
}



