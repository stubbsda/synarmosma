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
    /// This method returns as a rational number the arithmetic or harmonic (according to the third argument's value) mean of two integers.
    inline friend Rational compute_mean(int,int,std::string = "ARITHMETIC");
    /// This method returns as a rational number the arithmetic or harmonic (according to the third argument's value) mean of an integer and an instance of the Rational class. 
    inline friend Rational compute_mean(int,const Rational&,std::string = "ARITHMETIC");
    /// This method returns as a rational number the arithmetic or harmonic (according to the third argument's value) mean of an integer and an instance of the Rational class. 
    inline friend Rational compute_mean(const Rational&,int,std::string = "ARITHMETIC");
    /// This method returns as a rational number the arithmetic or harmonic (according to the third argument's value) mean of two instances of the Rational class. 
    inline friend Rational compute_mean(const Rational&,const Rational&,std::string = "ARITHMETIC");
    /// This method returns as a rational number the arithmetic or harmonic (according to the third argument's value) mean of a vector of integers.
    inline friend Rational compute_mean(const std::vector<int>&,std::string = "ARITHMETIC");
    /// This method returns as a rational number the arithmetic or harmonic (according to the third argument's value) mean of a vector of instances of the Rational class.
    inline friend Rational compute_mean(const std::vector<Rational>&,std::string = "ARITHMETIC");
    /// This method implements the unary negation operator for an instance of the Rational class.
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
    /// This method overrides the ostream operator so as to do a pretty print of an instance of the Rational class.
    friend std::ostream& operator <<(std::ostream&,const Rational&);
  };

  Rational compute_mean(int a,int b,std::string type)
  {
    return compute_mean(Rational(a),Rational(b),type);
  }

  Rational compute_mean(int a,const Rational& b,std::string type)
  {
    return compute_mean(Rational(a),b,type);
  }

  Rational compute_mean(const Rational& a,int b,std::string type)
  {
    return compute_mean(a,Rational(b),type);
  }

  Rational compute_mean(const Rational& a,const Rational& b,std::string type)
  {
    boost::to_upper(type);
    if (type != "ARITHMETIC" && type != "HARMONIC") throw std::invalid_argument("The type of the mean must be arithmetic or harmonic!");
    Rational output;
    if (type == "ARITHMETIC") {
      output = Rational(1,2)*(a + b);
    }
    else {
      // The harmonic mean...
      Rational q1 = a,q2 = b;
      q1.invert(); q2.invert();
      output = Rational(2,1)/(q1 + q2);
    }
    output.normalize();
    return output;
  }

  Rational compute_mean(const std::vector<int>& vx,std::string type)
  {
    unsigned int i,n = vx.size();
    std::vector<Rational> qx;

    for(i=0; i<n; ++i) {
      qx.push_back(Rational(vx[i]));
    }
    return compute_mean(qx,type);
  }

  Rational compute_mean(const std::vector<Rational>& vx,std::string type)
  {
    unsigned int n = vx.size();
    if (n < 2) throw std::invalid_argument("Computing the mean requires at least two numbers!");
    boost::to_upper(type);
    if (type != "ARITHMETIC" && type != "HARMONIC") throw std::invalid_argument("The type of the mean must be arithmetic or harmonic!");
    unsigned int i;
    Rational output(0);

    if (type == "ARITHMETIC") {
      for(i=0; i<n; ++i) {
        output = output + vx[i];
      }
      output = Rational(1,n)*output;
    }
    else {
      Rational q;
      for(i=0; i<n; ++i) {
        q = vx[i]; q.invert();
        output = output + q;
      }
      output = Rational(n,1)/output;
    }
    output.normalize();
    return output;
  }

  int convert(const Rational&,int);
}
#endif
