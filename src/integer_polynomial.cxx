#include "integer_polynomial.h"

using namespace SYNARMOSMA;

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

template<class kind>
Integer_Polynomial<kind>::Integer_Polynomial(const std::vector<kind>& t)
{
  if (t.empty()) return;
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
  if (t.empty()) return;

  degree = t.size() - 1;
  for(i=0; i<t.size(); ++i) {
    terms.push_back(t[i]);
  }

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
  /// This method is an instantiation of write_terms() for the case of a polynomial over multiprecision integers, needed so that the write_ZZ routine can be called to handle the possibility of very large numbers.
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
  /// This method is an instantiation of read_terms() for the case of a polynomial over multiprecision integers, needed so that the read_ZZ routine can be called to handle the possibility of very large numbers.
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

namespace SYNARMOSMA {
  template<>
  /// This method is an instantiation of the evaluate() method for the case of a multiprecision integer, needed to handle the modulus calculation with a multiprecision representation of the characteristic
  void Integer_Polynomial<NTL::ZZ>::scale_coefficients()
  {
    if (characteristic == 0) return;

    unsigned int i;
    NTL::ZZ q,base = NTL::to_ZZ(characteristic);
    const NTL::ZZ z(0);

    for(i=0; i<=degree; ++i) {
      q = terms[i];
      if (q < z) {
        do {
          q += base;
          if (q >= z) break;
        } while(true);
      }
      terms[i] = q % base;
    }
  }
}

template<class kind>
void Integer_Polynomial<kind>::scale_coefficients()
{
  if (characteristic == 0) return;

  unsigned int i;
  kind q;

  for(i=0; i<=degree; ++i) {
    q = terms[i];
    if (q < Integer_Polynomial<kind>::zero) {
      do {
        q += characteristic;
        if (q >= Integer_Polynomial<kind>::zero) break;
      } while(true);
    }
    terms[i] = q % characteristic;
  }
}

template<class kind>
void Integer_Polynomial<kind>::compute_irreducibility(long limit)
{
  // We will try to use Eisenstein's criterion on this polynomial...
  unsigned int i;
  long p;
  bool good,success = false;
  NTL::PrimeSeq s;
  
  do {
    p = s.next();
    if (p > limit) break;
    if (terms[degree] % p == 0) continue;
    good = true;
    for(i=0; i<degree; ++i) {
      if (terms[i] % p != 0) {
        good = false;
        break;
      }
    }
    if (!good) continue;
    if (terms[0] % (p*p) != 0) success = true;
  } while(!success);
  if (success) irreducible = true;
}

template<class kind>
bool Integer_Polynomial<kind>::compute_galois_group(Group* G) const
{
  G->clear();

  if (degree < 2) {
    G->initialize("ORDER",1);
    return true;
  }
  if (degree == 2) {
    G->initialize("SYMMETRIC",2);
    return true;
  }
  if (degree == 3) {
    if (!irreducible) return false;
    Rational s1(-terms[2],terms[3]),s2(terms[1],terms[3]),s3(-terms[0],terms[3]);
    Rational delta = s1*s1*s2*s2 + 18*s1*s2*s3 - 27*s3*s3 - 4*s1*s1*s1*s3 - 4*s2*s2*s2;
    // Now check if the numerator and denominator of delta are both perfect squares...
    if (delta.perfect_square()) {
      G->initialize("ALTERNATING",3);
    }
    else {
      G->initialize("SYMMETRIC",3);
    }
    return true;
  }

  if (irreducible && NTL::ProbPrime(degree)) {
    // If this polynomial has precisely two non-real roots, the Galois group is S(degree); 
    // however the degree of this polynomial must be at least five and using Descartes' rule 
    // of signs we can at best conclusively demonstrate the existence of only two real roots, 
    // one positive and one negative. 
    unsigned int i,ncomplex,nchange = 0,nreal = 0; 
    int s,nsign;

    nsign = (terms[0] > 0) ? 1 : -1;
    // First the positive real roots...
    for(i=1; i<=degree; ++i) {
      if (nsign*terms[i] < 0) {
        nchange++;
        nsign *= -1;
      }
    }
    if (nchange < 2) {
      nreal += nchange;
    }
    else if (nchange == 2) {
      kind y1,y2,x = zero;
      bool rfound = false;

      y1 = evaluate(x);
      for(i=1; i<=50; ++i) {
        x += unity;
        y2 = evaluate(x);
        if (y1*y2 < 0) {
          rfound = true;
          break;
        }
      }
      // If we find a positive root, that means in fact there are 
      // two of them...
      if (rfound) nreal += 2;
    }
    // Now the negative real roots...
    nchange = 0;
    nsign = (terms[0] > 0) ? 1 : -1;
    for(i=1; i<=degree; ++i) {
      s = ((i % 2) == 0) ? 1 : -1;
      if (nsign*s*terms[i] < 0) {
        nchange++;
        nsign *= -1;
      } 
    }
    if (nchange < 2) {
      nreal += nchange;
    }
    else if (nchange == 2) {
      kind y1,y2,x = zero;
      bool rfound = false;

      y1 = evaluate(x);
      for(i=1; i<=50; ++i) {
        x = x - unity;
        y2 = evaluate(x);
        if (y1*y2 < 0) {
          rfound = true;
          break;
        }
      }
      // If we find a negative root, that means in fact there are 
      // two of them...
      if (rfound) nreal += 2;
    }
    // Finally, the number of non-real roots is the degree minus the number of real roots...
    ncomplex = degree - nreal;
    // Except that the number of complex roots must be even since the polynomial is over the 
    // reals...
    if ((ncomplex % 2) != 0) ncomplex--;
    if (ncomplex == 2) {
      G->initialize("SYMMETRIC",degree);
      return true;
    }
  }

  return false;
}

template<class kind>
void Integer_Polynomial<kind>::simplify()
{
  unsigned int i,d = 0;

  scale_coefficients();
  for(i=0; i<=degree; ++i) {
    if (terms[i] == Integer_Polynomial<kind>::zero) continue;
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
  if (terms[degree] == Integer_Polynomial<kind>::unity) monic = true;
  if (terms[0] == Integer_Polynomial<kind>::zero) homogeneous = true;
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
kind Integer_Polynomial<kind>::get_value(unsigned int n) const
{
  if (n > degree) throw std::invalid_argument("The integer polynomial coefficient cannot exceed the degree!");
 
  return terms[n];
}

template<class kind>
void Integer_Polynomial<kind>::get_value(std::vector<kind>& vx) const
{
  vx = terms;
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

template<class kind>
void Integer_Polynomial<kind>::set_value(const std::vector<kind>& vx)
{
  degree = vx.size() - 1;
  terms = vx;
  simplify();
}

namespace SYNARMOSMA {
  template<>
  /// This method is an instantiation of the evaluate() method for the case of a multiprecision integer, needed to handle the modulus calculation with a multiprecision representation of the characteristic.
  NTL::ZZ Integer_Polynomial<NTL::ZZ>::evaluate(NTL::ZZ x) const
  {
    if (characteristic > 0) {
      if (x >= characteristic) throw std::invalid_argument("Argument exceeds field characteristic!");  
    }
    int i;
    NTL::ZZ y(0);

    for(i=degree; i>=0; --i) {
      y = y*x + terms[i];
    }
    if (characteristic > 0) y = y % NTL::to_ZZ(characteristic);    
    return y;
  }
}

template<class kind>
kind Integer_Polynomial<kind>::evaluate(kind x) const
{
  if (characteristic > 0) {
    if (x >= (signed) characteristic) throw std::invalid_argument("Argument exceeds field characteristic!");
  }
  int i;
  kind y = Integer_Polynomial<kind>::zero;
  // Use Horner's method to speed evaluation of the polynomial
  for(i=degree; i>0; --i) {
    y = y*x + terms[i];
  }
  if (characteristic > 0) y = y % characteristic;
  return y;
}

namespace SYNARMOSMA {
  template<>
  /// This method is an instantiation of the derivative() method, needed to handle the conversion of the formal differentiation expression for multiprecision integers.
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
  /// This method is an instantiation of the reduce() method, needed to handle the conversion of the NTL::ZZ coefficients to integers of "long" type.
  Integer_Polynomial<int> Integer_Polynomial<NTL::ZZ>::reduce(unsigned int p)
  {
    if (!NTL::ProbPrime(p)) throw std::invalid_argument("Integer_Polynomial must be reduced over a prime characteristic!");

    unsigned int i;
    long temp;
    NTL::ZZ q,base = NTL::to_ZZ(p);
    std::vector<int> nterms;
    const NTL::ZZ z(0);

    for(i=0; i<=degree; ++i) {
      q = terms[i];
      if (q < z) {
        do {
          q += base;
          if (q >= z) break;
        } while(true);
      }
      q = q % base;
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
  kind q;
  std::vector<int> nterms;

  for(i=0; i<=degree; ++i) {
    q = terms[i];
    if (q < Integer_Polynomial<kind>::zero) {
      do {
        q += p;
        if (q >= Integer_Polynomial<kind>::zero) break;
      } while(true);
    }
    q = q % p;
    nterms.push_back(q);
  }
  Integer_Polynomial<int> output(nterms,p);
  return output;
}

