#include "polynomial.h"

using namespace SYNARMOSMA;

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

namespace SYNARMOSMA {
  template<>
  /// This method is an instantiation of simplify() for the case of a polynomial over the rationals, needed due to the special representation of the number 1 (= 1/1) in the Rational class.
  void Polynomial<Rational>::simplify()
  {
    unsigned int i,d = 0;
    const Rational unity(1);

    for(i=0; i<=degree; ++i) {
      if (terms[i].is_null()) continue;
      d = i;
    }
    if (d < degree) {
      std::vector<Rational> nterms;

      for(i=0; i<=d; ++i) {
        nterms.push_back(terms[i]);
      }
      terms = nterms;
      degree = d;
    }
    if (terms[degree] == unity) monic = true;
    if (terms[0].is_null()) homogeneous = true;
  }
}

template<class kind>
void Polynomial<kind>::simplify()
{
  unsigned int i,d = 0;

  for(i=0; i<=degree; ++i) {
    if (std::abs(terms[i]) < std::numeric_limits<double>::epsilon()) continue;
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
  if (std::abs(terms[degree] - 1.0) == std::numeric_limits<double>::epsilon()) monic = true;
  if (std::abs(terms[0]) < std::numeric_limits<double>::epsilon()) homogeneous = true;
}

namespace SYNARMOSMA {
  template<>
  /// This method is an instantiation of write_terms() for the case of a polynomial over the rationals, needed so that the write_ZZ routine can be called to write the numerator and denominator.
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
  /// This method is an instantiation of read_terms() for the case of a polynomial over the rationals, needed so that the read_ZZ routine can be called to read the numerator and denominator.
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
void Polynomial<kind>::initialize(int L)
{
  if (L < 1) throw std::invalid_argument("The argument for Polynomial::initialize must be positive!");
  unsigned int i;
  kind test;
  Random RND;

  for(i=0; i<degree; ++i) {
    terms.push_back(kind(RND.irandom(-L,L)));
  }

  test = kind(RND.irandom(1,L));
  if (RND.irandom(2) == 0) test = -test;
  terms.push_back(test);

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
kind Polynomial<kind>::get_value(unsigned int n) const
{
  if (n > degree) throw std::invalid_argument("The polynomial coefficient cannot exceed the degree!");

  return terms[n];
}

template<class kind>
void Polynomial<kind>::get_value(std::vector<kind>& vx) const
{
  vx = terms;
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
kind Polynomial<kind>::evaluate(kind x) const
{
  kind y = kind(0);
  // Use Horner's method to speed evaluation of the polynomial
  for(int i=degree; i>=0; --i) {
    y = y*x + terms[i];
  }
  return y;
}

namespace SYNARMOSMA {
  template<>
  /// This method is an instantiation of derivative() for the case of a polynomial over the rationals, needed so that the coefficient of each term is correctly identified as a Rational instance.
  Polynomial<Rational> Polynomial<Rational>::derivative() const
  {
    unsigned int i;
    Rational q;
    Polynomial<Rational> output(degree-1);
    for(i=0; i<degree-1; ++i) {
      q = Rational(1+i,1);
      output.set_value(q*terms[i+1],i);
    }
    output.simplify();
    return output;
  }
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



