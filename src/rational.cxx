#include "rational.h"

using namespace SYNARMOSMA;

Rational::Rational()
{

}

Rational::Rational(int x)
{
  n = NTL::to_ZZ(x);
  compute_height();
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

Rational::Rational(double t,double tolerance)
{
  if (tolerance < std::numeric_limits<double>::epsilon()) throw std::invalid_argument("The tolerance for the rational approximation must be greater than zero!");
  int q;
  long double f,g,wcopy = std::abs(t);
  std::vector<int> cfractions;
  std::vector<int>::reverse_iterator rit;
  Rational sum;

  do {
    f = std::modf(wcopy,&g);
    q = int(g);
    wcopy = 1.0/f;
    // Now compute the rational...
    sum = Rational(q,1);
    for(rit=cfractions.rbegin(); rit!=cfractions.rend(); ++rit) {
      sum = 1/sum + Rational(*rit,1);
    }
    cfractions.push_back(q);
    if (std::abs(sum.to_double() - std::abs(t)) < tolerance) break;
  } while(true);

  n = sum.get_numerator();
  d = sum.get_denominator();
  if (t < -std::numeric_limits<double>::epsilon()) n = -n;    
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
  compute_height();
}

void Rational::normalize()
{
  // This method will normalize the sign of n and d, make n and d co-prime and finally compute the height.
  if (d < 0) {
    n = -n;
    d = -d;
  }
  NTL::ZZ c = NTL::GCD(NTL::abs(n),NTL::abs(d));
  if (c > 1) {
    n = n/c;
    d = d/c;
  }
  compute_height();
}

int Rational::serialize(std::ofstream& s) const
{
  int count = 0;

  count += write_ZZ(s,n);
  count += write_ZZ(s,d);

  return count;
}

int Rational::deserialize(std::ifstream& s)
{
  int count = 0;

  count += read_ZZ(s,n);
  count += read_ZZ(s,d);

  height = get_height();

  return count;
}

unsigned int Rational::reduce(unsigned int p) const
{
  if (d % p == 0) throw std::invalid_argument("The denominator must not be a multiple of the characteristic!");
  // This method must convert this rational number to an element of GF(p) for prime p, by first reducing 
  // the numerator and denominator over p and then doing the multiplication in the finite field.
  unsigned int np,dp;

  np = n % p;
  dp = d % p;
  if (np == 0 || dp == 1) return np;
  // In this more complex case, we need to actually carry out the multiplication np * dp^(-1) in GF(p); we 
  // begin by using the Euclidean algorithm to find the value of dp^(-1) in GF(p), i.e. solving the congruence
  // a*dp = 1 mod p for a
  unsigned int i = 1;
  std::vector<unsigned int> r,q;
  std::vector<int> s,t;
  s.push_back(1); s.push_back(0);
  t.push_back(0); t.push_back(1);
  r.push_back(p); r.push_back(dp);
  q.push_back(0);
  
  do {
    q.push_back(r[i-1]/r[i]);
    r.push_back(r[i-1] - q[i]*r[i]);
    s.push_back(s[i-1] - q[i]*s[i]);
    t.push_back(t[i-1] - q[i]*t[i]);
    if (r[i+1] == 1) break;
    i++;
  } while(true);
  // Correct the result so that it belongs to the right domain for this Galois field...
  int m,tf = t[t.size()-1];
  do {
    m = tf + p;
    if (m >= 0) break;
  } while(true);
  unsigned int um = (unsigned) m;
  // Now do the multiplication and return the result...
  unsigned int output = (np * um) % p;
  
  return output; 
}

bool Rational::perfect_square() const
{
  // We use binary search to check if both the numerator and 
  // denominator are perfect squares, returning true if so. 
  bool output = true;
  NTL::ZZ lo,hi,mid,sq,nt;

  nt = (n < 0) ? -n : n;
  if (nt > 1) {
    output = false; 
    lo = NTL::to_ZZ(2);
    hi = NTL::to_ZZ(nt);
    do {
      mid = lo + (hi - lo)/2;
      sq = mid*mid;
      if (sq == nt) {
        output = true;
        break;
      }
      if (sq < nt) {
        lo = mid;
      }
      else {
        hi = mid;
      }
    } while((hi - lo) > 1);
    if (!output) return false;
  }

  nt = (d < 0) ? -d : d;
  if (nt > 1) {
    output = false;
    lo = NTL::to_ZZ(2);
    hi = NTL::to_ZZ(nt);
    do {
      mid = lo + (hi - lo)/2;
      sq = mid*mid;
      if (sq == nt) {
        output = true;
        break;
      }
      if (sq < nt) {
        lo = mid;
      }
      else {
        hi = mid;
      }
    } while((hi - lo) > 1);
  }
  return output;
}

long Rational::agreeableness() const
{
  // First calculate the least common multiple M of the numerator 
  // and denominator, then compute the prime decomposition of M. 
  // This method then returns 
  // sigma - nf + 1
  // where sigma is the sum of all the factors and nf the number 
  // of factors, counted with multiplicity. For example, given the 
  // ratio 18/7 so that M = 18*7 = 2*3*3*7, the agreeableness will 
  // be (2 + 3 + 3 + 7) - 4 + 1 = 12. According to Euler, the lower 
  // the number, the more agreeable the pitch ratio.
  if (n == 0) return -1;

  long i,nf,sigma = 0;
  NTL::ZZ c = NTL::GCD(NTL::abs(n),NTL::abs(d));
  std::vector<NTL::ZZ> factors;
  NTL::PrimeSeq s;
  NTL::ZZ L,q,p = (n*d)/c;

  if (p < 0) p = -p;

  if (NTL::ProbPrime(p)) {
    factors.push_back(p);
  }
  else {
    L = NTL::to_ZZ(long(std::sqrt(NTL::to_long(p))));
    do {
      q = NTL::to_ZZ(s.next());
      do {
        if ((p % q) != 0) break;
        factors.push_back(q);
        p = p/q;
      } while(true);
      if (q > L) break;
    } while(true);
  }

  nf = (signed) factors.size();
  for(i=0; i<nf; ++i) {
    sigma += NTL::to_long(factors[i]);
  }

  return sigma - nf + 1;
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

  bool operator >=(const Rational& r1,int alpha)
  {
    Rational r2(alpha);
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

  bool operator <=(const Rational& r1,int alpha)
  {
    Rational r2(alpha);
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

    r = q.get_numerator() % p;
    if (q.get_denominator() == 1) {
      for(i=0; i<p; ++i) {
        if (r == i) {
          output = i;
          break;
        }
      }
      return output;
    }
    s = q.get_denominator() % p;
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
