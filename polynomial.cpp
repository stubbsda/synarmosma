#include "polynomial.h"

using namespace SYNARMOSMA;

extern Random RND;

template<class kind>
Polynomial<kind>::Polynomial()
{
  characteristic = 0;
  degree = RND.irandom(5,10);
  initialize();
}

template<class kind>
Polynomial<kind>::Polynomial(unsigned int n)
{
  characteristic = 0;
  degree = n;
  initialize();
}

template<class kind>
Polynomial<kind>::Polynomial(unsigned int n,unsigned int p)
{
  characteristic = p;
  degree = n;
  initialize();
}

template<class kind>
Polynomial<kind>::Polynomial(const std::vector<kind>& t)
{
  unsigned int i;

  degree = t.size()-1;
  for(i=0; i<t.size(); ++i) {
    terms.push_back(t[degree-i]);
  }
  normed = false;
  if (terms[degree] == kind(1)) normed = true;
}

template<class kind>
Polynomial<kind>::Polynomial(const Polynomial& source)
{
  unsigned int i;

  terms.clear();
  degree = source.degree;
  characteristic = source.characteristic;
  normed = source.normed;
  homogeneous = source.homogeneous;
  irreducible = source.irreducible;
  for(i=0; i<=degree; ++i) {
    terms.push_back(source.terms[i]);
  }
}

template<class kind>
Polynomial<kind>::~Polynomial()
{

}

template<class kind>
void Polynomial<kind>::initialize()
{
  unsigned int i;
  int test;
  bool flag = true;

  if (characteristic == 0) {
    for(i=0; i<degree; ++i) {
      terms.push_back(RND.irandom(-25,25));
    }
    do {
      test = RND.irandom(-25,25);
      if (test != 0) flag = false;
    } while(flag);
    terms.push_back(test);
  }
  else {
    for(i=0; i<degree; ++i) {
      terms.push_back(RND.irandom(characteristic));
    }
    terms.push_back(RND.irandom(1,characteristic-1));
  }
  normed = false;
  if (terms[degree] == kind(1)) normed = true;
  homogeneous = false;
  if (terms[0] == kind(0)) homogeneous = true;
  irreducible = true;
}

template<class kind>
Polynomial<kind>& Polynomial<kind>::operator =(const Polynomial<kind>& source)
{
  if (this == &source) return *this;
  terms.clear();
  unsigned int i;
  degree = source.degree;
  characteristic = source.characteristic;
  normed = source.normed;
  homogeneous = source.homogeneous;
  irreducible = source.irreducible;
  for(i=0; i<=degree; ++i) {
    terms.push_back(source.terms[i]);
  }
  return *this;
}

template<class kind>
kind Polynomial<kind>::get_value(unsigned int i) const
{
  return terms[i];
}

template<class kind>
void Polynomial<kind>::set_value(kind x,unsigned int i)
{
  if (i <= degree) {
    terms[i] = x;
    return;
  }
}

namespace SYNARMOSMA {
  template<>
  unsigned int Polynomial<unsigned int>::evaluate(unsigned int x)
  {
    if (x > characteristic) {
      std::cerr << "Cannot evaluate this argument!" << std::endl;
      std::exit(1);
    }
    int i;
    unsigned int y = 0;
    for(i=degree; i>=0; --i) {
      y = y*x + terms[i];
    }
    y = y % characteristic;
    return y;
  }
}

template<class kind>
kind Polynomial<kind>::evaluate(kind x)
{
  int i;
  kind y = 0;
  // Use Horner's method to speed evaluation of the polynomial
  for(i=degree; i>=0; --i) {
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
  return output;
}

namespace SYNARMOSMA {
  template <>
  Polynomial<unsigned int> Polynomial<std::complex<double> >::reduce(unsigned int p)
  {
    Polynomial<unsigned int> output(degree);
    return output;
  }

  template<>
  Polynomial<unsigned int> Polynomial<double>::reduce(unsigned int p)
  {
    Polynomial<unsigned int> output(degree);
    return output;
  }
}

template<class kind>
Polynomial<unsigned int> Polynomial<kind>::reduce(unsigned int p)
{
  std::vector<unsigned int> new_terms;
  unsigned int i;

  for(i=0; i<=degree; ++i) {
    new_terms.push_back(convert(terms[i],p));
  }
  Polynomial<unsigned int> output(new_terms);
  return output;
}

template<class kind>
void Polynomial<kind>::factorize(std::vector<Polynomial>& output)
{

}

template<class kind>
Polynomial<kind> operator -(const Polynomial<kind>& p1,const Polynomial<kind>& p2)
{
  std::vector<kind> new_terms;
  unsigned int i,mu = minimum(p1.degree,p2.degree);
  for(i=0; i<mu; ++i) {
    new_terms.push_back(p1.terms[i] - p2.terms[i]);
  }
  if (p1.degree > p2.degree) {
    for(i=mu; i<p1.degree; ++i) {
      new_terms.push_back(p1.terms[i]);
    }
  }
  else if (p2.degree > p1.degree) {
    for(i=mu; i<p2.degree; ++i) {
      new_terms.push_back(-p2.terms[i]);
    }
  }
  Polynomial<kind> output(new_terms);
  return output;
}

template<class kind>
Polynomial<kind> operator +(const Polynomial<kind>& p1,const Polynomial<kind>& p2)
{
  std::vector<kind> new_terms;
  unsigned int i,mu = minimum(p1.degree,p2.degree);
  for(i=0; i<mu; ++i) {
    new_terms.push_back(p1.terms[i] + p2.terms[i]);
  }
  if (p1.degree > p2.degree) {
    for(i=mu; i<p1.degree; ++i) {
      new_terms.push_back(p1.terms[i]);
    }
  }
  else if (p2.degree > p1.degree) {
    for(i=mu; i<p2.degree; ++i) {
      new_terms.push_back(p2.terms[i]);
    }
  }
  Polynomial<kind> output(new_terms);
  return output;
}


