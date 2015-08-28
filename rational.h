/*
  Copyright 2014 Daniel Stubbs

  This file is part of Synarmosma.

  Synarmosma is free software: you can redistribute it and/or modify 
  it under the terms of the GNU General Public License as published by 
  the Free Software Foundation, either version 3 of the License, or 
  (at your option) any later version.

  Synarmosma is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with Synarmosma.  If not, see <http://www.gnu.org/licenses/>.
*/

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
