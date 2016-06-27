#include "global.h"

#ifndef _rationalh
#define _rationalh

// A class for rational numbers, i.e. a number of the form n/d where 
// n and d are whole numbers. We assume that d is always greater 
// than zero.
namespace SYNARMOSMA {  
  class Rational {
   private:
    NTL::ZZ d,n;
    double height;
  
    void normalize();
    void invert();
   public:
    Rational();
    Rational(signed int); 
    Rational(signed int,signed int);
    Rational(const NTL::ZZ&,const NTL::ZZ&);
    Rational& operator =(const Rational&);
    Rational operator -();
    Rational(const Rational&);
    ~Rational();
    NTL::ZZ numerator() const;
    NTL::ZZ denominator() const;
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

  unsigned int convert(const Rational&,unsigned int);
  Rational qdiv(const Rational&,const Rational&);
}
#endif 
