#include "global.h"

#ifndef _rationalh
#define _rationalh

namespace SYNARMOSMA {  
  /// A class for rational numbers, i.e. a number of the form n/d where n and d are whole numbers; we assume that d is always greater than zero.
  class Rational {
   private:
    /// The numerator of a rational number, stored as a multi-precision integer using the NTL::ZZ type.
    NTL::ZZ n;
    /// The denominator of a rational number, stored as a multi-precision integer using the NTL::ZZ type.
    NTL::ZZ d;
    /// The height of this rational number, defined to be log(abs(n)) when n >= d, otherwise log(abs(d)).
    double height;
  
    /// This method makes the numerator and denominator coprime and computes the height property of this instance.
    void normalize();
    /// Turns the rational q into its reciprocal 1/q, i.e. n => d', d => n'.
    void invert();
   public:
    Rational();
    Rational(int); 
    Rational(int,int);
    Rational(const NTL::ZZ&,const NTL::ZZ&);
    Rational& operator =(const Rational&);
    Rational operator -();
    Rational(const Rational&);
    ~Rational();
    inline NTL::ZZ numerator() const {return n;};
    inline NTL::ZZ denominator() const {return d;};
    friend Rational operator -(const Rational&);
    friend Rational operator +(const Rational&,const Rational&);
    friend Rational operator -(const Rational&,const Rational&);
    friend Rational operator *(const Rational&,const Rational&);
    friend Rational operator /(const Rational&,const Rational&); 
    friend bool operator ==(const Rational&,const Rational&);
    friend bool operator ==(const Rational&,int);
    friend bool operator !=(const Rational&,const Rational&);
    friend bool operator !=(const Rational&,int);
    friend bool operator <=(const Rational&,const Rational&);
    friend bool operator >=(const Rational&,const Rational&);
    friend bool operator <(const Rational&,const Rational&);
    friend bool operator <(const Rational&,int);
    friend bool operator >(const Rational&,const Rational&);
    friend bool operator >(const Rational&,int);
    friend std::ostream& operator <<(std::ostream&,const Rational&);
  };

  int convert(const Rational&,int);
}
#endif 
