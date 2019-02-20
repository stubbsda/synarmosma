#include "rational.h"

using namespace SYNARMOSMA;

Rational::Rational()
{
  n = 0;
  d = 1;
}

Rational::Rational(int x)
{
  n = NTL::to_ZZ(x);
  d = 1;
}

Rational::Rational(int x,int y)
{
  if (y == 0) throw std::invalid_argument("The denominator of a rational cannot be zero!");
  n = NTL::to_ZZ(x);
  d = NTL::to_ZZ(y);
  normalize();
}

Rational::Rational(const NTL::ZZ& x,const NTL::ZZ& y)
{
  if (y == 0) throw std::invalid_argument("The denominator of a rational cannot be zero!");
  n = x;
  d = y;
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

void Rational::invert()
{
  if (n == 0) throw std::runtime_error("Cannot invert a rational equal to zero!");
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
  height = (NTL::abs(n) > NTL::abs(d)) ? NTL::log(NTL::abs(n)) : NTL::log(NTL::abs(d));
}

namespace SYNARMOSMA {

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
    if (r2.n == 0) throw std::invalid_argument("Cannot divide by a rational equal to zero!");
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

  Rational operator -(const Rational& r)
  {
    Rational output = r; 
    if (output.d < 0) output.d = -output.d;
    if (output.n < 0) output.n = -output.n;
    return output;
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

  bool operator >(const Rational& r1,int alpha)
  {
    Rational r2(alpha);
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

  bool operator <(const Rational& r1,int alpha)
  {
    Rational r2(alpha);
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

  bool operator !=(const Rational& r1,int alpha)
  {
    Rational r2(alpha);
    if (r1.n != r2.n || r1.d != r2.d) return true;
    return false;
  }

  bool operator ==(const Rational& r1,const Rational& r2)
  {
    if (r1.n == r2.n && r1.d == r2.d) return true;
    return false;
  }

  bool operator ==(const Rational& r1,int alpha)
  {
    Rational r2(alpha);
    if (r1.n == r2.n && r1.d == r2.d) return true;
    return false;
  }

  int convert(const Rational& q,int p)
  {
    // Here we have to use the relation q = r/s => s*q = r, i.e. 
    // we would like to find out which element q of GF(p) satisfies 
    // this equation.
    // First we will convert r and s from Z to Z/p...
    int i,output = 0;
    NTL::ZZ r,s,in1;

    if (!NTL::ProbPrime(p)) throw std::invalid_argument("Field characteristic must be prime!");

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
