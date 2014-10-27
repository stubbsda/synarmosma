#include "rational.h"

using namespace SYNARMOSMA;

Rational::Rational()
{
  n = 0;
  d = 1;
}

Rational::Rational(signed int x)
{
  n = NTL::to_ZZ(x);
  d = 1;
}

Rational::Rational(signed int x,signed int y)
{
  n = NTL::to_ZZ(x);
  d = NTL::to_ZZ(y);
  normalize();
}

Rational::Rational(const Rational& source)
{
  n = source.n;
  d = source.d;
  height = source.height;
}

Rational& Rational::operator =(const Rational& source)
{
  if (this == &source) return *this;
  n = source.n;
  d = source.d;
  height = source.height;
  return *this;
}

Rational Rational::operator -()
{
  Rational output;
  output.n = -n;
  output.d = d;
  output.height = height;
  return output;
}

Rational::~Rational()
{

}

NTL::ZZ Rational::numerator() const
{
  return n;
}

NTL::ZZ Rational::denominator() const
{
  return d;
}

void Rational::invert()
{
  NTL::ZZ temp = n;
  n = d;
  d = temp;
}

void Rational::normalize()
{
  // This method will make n and d coprime.
  NTL::ZZ c;
  c = NTL::GCD(n,d);
  if (c != 1) {
    n = n/c;
    d = d/c;
  }
  height = (abs(n) > abs(d)) ? log(n) : log(d);
}

namespace SYNARMOSMA {
  Rational abs(Rational source)
  {
    Rational output;
    output.n = abs(source.n);
    output.d = abs(source.d);
    return output;
  }

  Rational operator +(const Rational& r1,const Rational& r2)
  {
    Rational output(0);
    if (r1.d == r2.d) {
      output.d = r1.d;
      output.n = r1.n + r2.n;
      output.normalize();
      return output;
    }
    output.d = r1.d*r2.d;
    output.n = (r1.n*r2.d + r2.n*r1.d);
    output.normalize();
    return output;
  }

  Rational operator -(const Rational& r1,const Rational& r2)
  {
    Rational output(0);
    if (r1.d == r2.d) {
      output.d = r1.d;
      output.n = r1.n - r2.n;
      output.normalize();
      return output;
    }
    output.d = r1.d*r2.d;
    output.n = (r1.n*r2.d - r2.n*r1.d);
    output.normalize();
    return output;
  }

  Rational operator *(const Rational& r1,const Rational& r2)
  {
    Rational output(0);
    output.n = r1.n*r2.n;
    output.d = r1.d*r2.d;
    output.normalize();
    return output;
  }

  Rational operator /(const Rational& r1,const Rational& r2)
  {
    Rational output(0);
    output.n = r1.n*r2.d;
    output.d = r1.d*r2.n;
    output.normalize();
    return output;
  }

  std::ostream& operator <<(std::ostream& s,const Rational& source)
  {
    s << "(" << source.n << ")/(" << source.d << ")";
    return s;
  }

  bool operator >=(const Rational& r1,const Rational& r2)
  {
    if (r1.d == r2.d) {
      return (r1.n >= r2.n);
    }
    return ((r1.n*r2.d) >= (r2.n*r1.d));
  }

  bool operator >(const Rational& r1,const Rational& r2)
  {
    if (r1.d == r2.d) {
      return (r1.n > r2.n);
    }
    return ((r1.n*r2.d) > (r2.n*r1.d));
  }

  bool operator <(const Rational& r1,const Rational& r2)
  {
    if (r1.d == r2.d) {
      return (r1.n < r2.n);
    }
    return ((r1.n*r2.d) < (r2.n*r1.d));
  }

  bool operator <=(const Rational& r1,const Rational& r2)
  {
    if (r1.d == r2.d) {
      return (r1.n <= r2.n);
    }
    return ((r1.n*r2.d) <= (r2.n*r1.d));
  }

  bool operator !=(const Rational& r1,const Rational& r2)
  {
    if (r1.n != r2.n || r1.d != r2.d) return true;
    return false;
  }

  bool operator ==(const Rational& r1,const Rational& r2)
  {
    if (r1.n == r2.n && r1.d == r2.d) return true;
    return false;
  }

  Rational qdiv(const Rational& a,const Rational& b)
  {
    return a;
  }

  unsigned int convert(const Rational& q,unsigned int p)
  {
    // Here we have to use the relation q = r/s => s*q = r, i.e. 
    // we would like to find out which element q of GF(p) satisfies 
    // this equation.
    // First we will convert r and s from Z to Z/p...
    unsigned int i,output = 0;
    NTL::ZZ r,s,in1;

    r = q.numerator() % p;
    if (q.denominator() == 1) {
      for(i=0; i<p; ++i) {
        if (r == i) {
          output = i;
          break;
        }
      }
      return output;
    }
    s = q.denominator() % p;
    // This method isn't very efficient but it works.
    for(i=0; i<p; ++i) {
      in1 = s*i;
      if ((in1 % p) == r) {
        output = i;
        break;
      }
    }
    return output;
  }
}
