#include "polynomial.h"

using namespace SYNARMOSMA;

extern Random RND;

template<class kind>
Polynomial<kind>::Polynomial()
{
  degree = RND.irandom(5,10);
  initialize();
}

template<class kind>
Polynomial<kind>::Polynomial(unsigned int n)
{
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

namespace SYNARMOSMA {
  template<>
  void Polynomial<NTL::ZZ>::property_check()
  {
    if (terms[degree] == NTL::to_ZZ(1)) normed = true;
    if (terms[0] == NTL::to_ZZ(0)) homogeneous = true;
  }
}

template<class kind>
void Polynomial<kind>::property_check()
{
  if (terms[degree] == kind(1)) normed = true;
  if (terms[0] == kind(0)) homogeneous = true;
}

template<class kind>
Polynomial<kind>::Polynomial(const std::vector<kind>& t)
{
  unsigned int i;

  degree = t.size() - 1;
  for(i=0; i<t.size(); ++i) {
    terms.push_back(t[i]);
  }

  property_check();
}

template<class kind>
Polynomial<kind>::Polynomial(const std::vector<kind>& t,unsigned int p)
{
  unsigned int i;

  characteristic = p;
  degree = t.size() - 1;
  for(i=0; i<t.size(); ++i) {
    terms.push_back(t[i]);
  }
  irreducible = false;
  homogeneous = false;
  normed = false;
  property_check();
}

template<class kind>
Polynomial<kind>::Polynomial(const Polynomial& source)
{
  degree = source.degree;
  characteristic = source.characteristic;
  normed = source.normed;
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
  characteristic = 0;
  degree = 0;
  homogeneous = false;
  irreducible = false;
  normed = false;
}

namespace SYNARMOSMA {
  template<>
  int Polynomial<Rational>::write_terms(std::ofstream& s) const
  {
    unsigned int i;
    int count = 0;

    for(i=0; i<=degree; ++i) {
      count += write_ZZ(s,terms[i].numerator());
      count += write_ZZ(s,terms[i].denominator());
    }
    return count;
  }

  template<>
  int Polynomial<NTL::ZZ>::write_terms(std::ofstream& s) const
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
  s.write((char*)(&characteristic),sizeof(int)); count += sizeof(int);
  s.write((char*)(&homogeneous),sizeof(bool)); count += sizeof(bool);
  s.write((char*)(&irreducible),sizeof(bool)); count += sizeof(bool);
  s.write((char*)(&normed),sizeof(bool)); count += sizeof(bool);
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

  template<>
  int Polynomial<NTL::ZZ>::read_terms(std::ifstream& s)
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
  s.read((char*)(&characteristic),sizeof(int)); count += sizeof(int);
  s.read((char*)(&homogeneous),sizeof(bool)); count += sizeof(bool);
  s.read((char*)(&irreducible),sizeof(bool)); count += sizeof(bool);
  s.read((char*)(&normed),sizeof(bool)); count += sizeof(bool);
  count += read_terms(s);

  return count;
}

namespace SYNARMOSMA {
  template<>
  void Polynomial<NTL::ZZ>::initialize()
  {
    unsigned int i;
    NTL::ZZ test;
    bool flag = true;
    const NTL::ZZ zero = NTL::to_ZZ(0);

    for(i=0; i<degree; ++i) {
      terms.push_back(NTL::to_ZZ(RND.irandom(-25,25)));
    }
    do {
      test = NTL::to_ZZ(RND.irandom(-25,25));
      if (test != zero) flag = false;
    } while(flag);
    terms.push_back(test);

    irreducible = true;
    property_check();
  }
}

template<class kind>
void Polynomial<kind>::initialize()
{
  unsigned int i;
  kind test;
  bool flag = true;
  const kind zero = kind(0);

  if (characteristic == 0) {
    for(i=0; i<degree; ++i) {
      terms.push_back(kind(RND.irandom(-25,25)));
    }
    do {
      test = kind(RND.irandom(-25,25));
      if (test != zero) flag = false;
    } while(flag);
    terms.push_back(test);
  }
  else {
    for(i=0; i<degree; ++i) {
      terms.push_back(kind(RND.irandom(characteristic)));
    }
    terms.push_back(RND.irandom(1,characteristic-1));
  }
  irreducible = true;
  property_check();
}

template<class kind>
Polynomial<kind>& Polynomial<kind>::operator =(const Polynomial<kind>& source)
{
  if (this == &source) return *this;

  degree = source.degree;
  characteristic = source.characteristic;
  normed = source.normed;
  homogeneous = source.homogeneous;
  irreducible = source.irreducible;
  terms = source.terms;

  return *this;
}

template<class kind>
kind Polynomial<kind>::get_value(unsigned int i) const
{
  return terms[i];
}

template<class kind>
bool Polynomial<kind>::set_value(kind x,unsigned int i)
{
  if (i <= degree) {
    terms[i] = x;
    property_check();
    return true;
  }
  return false;
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

  template<>
  NTL::ZZ Polynomial<NTL::ZZ>::evaluate(NTL::ZZ x)
  {
    int i;
    NTL::ZZ y = NTL::to_ZZ(0);
    for(i=degree; i>=0; --i) {
      y = y*x + terms[i];
    }
    return y;
  }
}

template<class kind>
kind Polynomial<kind>::evaluate(kind x)
{
  int i;
  kind y = 0;
  // Use Horner's method to speed evaluation of the polynomial
  for(i=degree; i>0; --i) {
    y = y*x + terms[i];
  }
  return y;
}

namespace SYNARMOSMA {
  template<>
  Polynomial<NTL::ZZ> Polynomial<NTL::ZZ>::derivative() const
  {
    unsigned int i;
    Polynomial<NTL::ZZ> output(degree-1);
    for(i=0; i<degree-1; ++i) {
      output.set_value(NTL::to_ZZ(1+i)*terms[i+1],i);
    }
    output.property_check();
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
  output.property_check();
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

  template<>
  Polynomial<unsigned int> Polynomial<NTL::ZZ>::reduce(unsigned int p)
  {
    unsigned int i;
    long temp;
    NTL::ZZ q,base = NTL::to_ZZ(p);
    std::vector<unsigned int> nterms;

    for(i=0; i<=degree; ++i) {
      q = terms[i] % base;
      NTL::conv(temp,q);
      nterms.push_back((unsigned int) temp);
    }
    Polynomial<unsigned int> output(nterms,p);
    return output;
  }
}

template<class kind>
Polynomial<unsigned int> Polynomial<kind>::reduce(unsigned int p)
{
  unsigned int i;
  std::vector<unsigned int> nterms;

  for(i=0; i<=degree; ++i) {
    nterms.push_back(convert(terms[i],p));
  }
  Polynomial<unsigned int> output(nterms,p);
  return output;
}



