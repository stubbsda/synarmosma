#include "global.h"

#ifndef _rationalh
#define _rationalh

namespace SYNARMOSMA {  
  /// A class representing rational numbers, i.e. a number of the form n/d where n and d are whole numbers. 
  class Rational {
   private:
    /// The numerator of a rational number q=n/d with n and d co-prime and d > 0, stored as a 
    /// multi-precision integer using the NTL::ZZ type.
    NTL::ZZ n = NTL::to_ZZ(0);
    /// The denominator of a rational number q=n/d with n and d co-prime and d > 0, stored as a 
    /// multi-precision integer using the NTL::ZZ type.
    NTL::ZZ d = NTL::to_ZZ(1);
    /// The height of the rational number q = n/d, defined to be log(abs(n)) when abs(q) > 1, 
    /// otherwise log(abs(d)).
    double height = 0.0;
  
    /// This method makes the numerator and denominator coprime and computes the height property of this instance.
    void normalize();
    /// This method turns the rational q into its reciprocal 1/q, i.e. n => d', d => n'.
    void invert();
    /// This method calculates the height property of the rational number. 
    inline void compute_height() {height = (NTL::abs(n) > NTL::abs(d)) ? NTL::log(NTL::abs(n)) : NTL::log(NTL::abs(d));};
   public:
    /// The default constructor which sets the properties to their default value and exits.
    Rational();
    /// This is the constructor for a rational integer, i.e. q = n/1 where n is the constructor's argument.
    Rational(int); 
    /// This constructor is for a standard rational number, with the first argument set to the numerator and the second to the denominator, after which it calls the normalize method.
    Rational(int,int);
    /// A constructor for a standard rational number but which directly accepts multi-precision integers for the numerator and denominator.
    Rational(const NTL::ZZ&,const NTL::ZZ&);
    /// The assignment operator for instances of the Rational class. 
    Rational& operator =(const Rational&);
    /// The unary negation operator which multiplies the numerator by -1.
    Rational operator -();
    /// The copy constructor for this class.
    Rational(const Rational&);
    /// A destructor that is empty.
    ~Rational();
    /// This method returns the "suavitas" (agreeableness) of a ratio of two integers, as a pitch ratio, if this notion makes sense in this case; otherwise the method returns -1. The more consonant a pitch ratio is musically, the lower the value returned by this method, with the optimal value of 1 for the ratio 1:1.
    long agreeableness() const;
    /// This method returns true if this rational number is zero and false otherwise.
    inline bool is_null() const {return (n == NTL::to_ZZ(0));};
    /// This method returns the numerator.
    inline NTL::ZZ get_numerator() const {return n;};
    /// This method returns the denominator.
    inline NTL::ZZ get_denominator() const {return d;};
    /// This method returns the height.
    inline double get_height() const {return height;};
    /// Another method to implement the unary negation operator for this class.
    friend Rational operator -(const Rational&);
    /// This method implements the addition operator for two instances of the Rational class.
    friend Rational operator +(const Rational&,const Rational&);
    /// This method implements the subtraction operator for two instances of the Rational class.
    friend Rational operator -(const Rational&,const Rational&);
    /// This method implements the multiplication operator for two instances of the Rational class.
    friend Rational operator *(const Rational&,const Rational&);
    /// This method implements the division operator for two instances of the Rational class.
    friend Rational operator /(const Rational&,const Rational&);
    /// One of a set of overloaded Boolean operators for the Rational class, testing order relations among rationals and integers. 
    friend bool operator ==(const Rational&,const Rational&);
    /// One of a set of overloaded Boolean operators for the Rational class, testing order relations among rationals and integers.
    friend bool operator ==(const Rational&,int);
    /// One of a set of overloaded Boolean operators for the Rational class, testing order relations among rationals and integers.
    friend bool operator !=(const Rational&,const Rational&);
    /// One of a set of overloaded Boolean operators for the Rational class, testing order relations among rationals and integers.
    friend bool operator !=(const Rational&,int);
    /// One of a set of overloaded Boolean operators for the Rational class, testing order relations among rationals and integers.
    friend bool operator <=(const Rational&,const Rational&);
    /// One of a set of overloaded Boolean operators for the Rational class, testing order relations among rationals and integers.
    friend bool operator <=(const Rational&,int);
    /// One of a set of overloaded Boolean operators for the Rational class, testing order relations among rationals and integers.
    friend bool operator >=(const Rational&,const Rational&);
    /// One of a set of overloaded Boolean operators for the Rational class, testing order relations among rationals and integers.
    friend bool operator >=(const Rational&,int);
    /// One of a set of overloaded Boolean operators for the Rational class, testing order relations among rationals and integers.
    friend bool operator <(const Rational&,const Rational&);
    /// One of a set of overloaded Boolean operators for the Rational class, testing order relations among rationals and integers.
    friend bool operator <(const Rational&,int);
    /// One of a set of overloaded Boolean operators for the Rational class, testing order relations among rationals and integers.
    friend bool operator >(const Rational&,const Rational&);
    /// One of a set of overloaded Boolean operators for the Rational class, testing order relations among rationals and integers.
    friend bool operator >(const Rational&,int);
    /// This method overrides the ostream operator so as to do a pretty print of an instance of the class.
    friend std::ostream& operator <<(std::ostream&,const Rational&);
  };

  int convert(const Rational&,int);
}
#endif 
