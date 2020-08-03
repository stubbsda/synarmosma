#include "random.h"
#include "rational.h"
#include "integer_polynomial.h"

#ifndef _polynomialh
#define _polynomialh

namespace SYNARMOSMA {
  template<class kind>
  class Polynomial;

  template<class kind>
  std::ostream& operator <<(std::ostream&,const Polynomial<kind>&);

  template<class kind>
  bool operator ==(const Polynomial<kind>&,const Polynomial<kind>&);

  template<class kind>
  Polynomial<kind> operator +(const Polynomial<kind>&,const Polynomial<kind>&);

  template<class kind>
  Polynomial<kind> operator -(const Polynomial<kind>&,const Polynomial<kind>&);

  template<class kind>
  Polynomial<kind> operator *(kind,const Polynomial<kind>&);

  template<class kind>
  Polynomial<kind> operator *(const Polynomial<kind>&,const Polynomial<kind>&);

  template<class kind>
  /// A class representing a low-degree polynomial over a floating point domain, modelling a polynomial over the real or complex numbers.
  class Polynomial {
   protected:
    /// This non-negative property is the polynomial's degree, and thus also 
    /// determines the length of Polynomial::terms which is equal to one plus 
    /// the degree.
    unsigned int degree = 0;
    /// This Boolean property is true if the polynomial cannot be factorized over its 
    /// domain, i.e. it has no solutions over this domain; if the base type is the complex 
    /// numbers then of course this is always false.
    bool irreducible = false;
    /// This Boolean property is true if the polynomial has no constant term, i.e. 
    /// \f$p(x) = a_d x^d + \dots + a_2 x^2 + a_1 x\f$ where \f$d\f$ is the degree.
    bool homogeneous = false;
    /// This Boolean property is true if the coefficient of the highest term in the polynomial 
    /// is one, so that \f$ p(x) = x^d + a_{d-1}x^{d-1} + \dots + a_1 x + a_0\f$ where \f$d\f$ 
    /// is the degree.
    bool monic = false;
    /// This is the principal property of the class and contains a list of the polynomial coefficients 
    /// stored in a dense manner, so that this class is intended for low-degree polynomials and not one 
    /// like \f$x^{250} - 7x +1\f$. The length of the vector is always Polynomial::degree plus one.
    std::vector<kind> terms;
  
    /// This method constructs the polynomial \f$p(x) = a_d x^d + a_{d-1}x^{d-1} + \dots + a_1 x + a_0\f$ where \f$d\f$ is Polynomial::degree and the \f$a_i\f$ (\f$a_d\ne 0\f$) are uniform random variates on the interval \f$[-L,L]\f$ where \f$L\f$ is the method's argument.
    void initialize(int = 25);
    /// This method verifies that the degree corresponds to the highest non-zero coefficient and sets the values of Polynomial::monic and Polynomial::homogeneous.
    void simplify();
    /// This method writes the content of the polynomial itself in binary format to an output stream and is used by the serialize() method; it returns the number of bytes written.
    int write_terms(std::ofstream&) const;
    /// This method reads the content of the polynomial itself in binary format from an input stream is used by the deserialize() method; it returns the number of bytes read.
    int read_terms(std::ifstream&);
   public:
    /// The default constructor which does nothing.
    Polynomial();
    /// This constructor accepts as its argument the degree \f$d\f$ of the polynomial and constructs the polynomial \f$x^d + x^{d-1} + \dots + x + 1\f$.
    Polynomial(unsigned int);
    /// This constructor uses its argument as the value for Polynomial::terms and sets the other class properties appropriately by calling the simplify() method.
    Polynomial(const std::vector<kind>&);
    /// The destructor which does nothing.
    ~Polynomial();
    /// The overloaded assignment operator for this class, which first calls clear() and then behaves exactly like the copy constructor for this class.
    Polynomial& operator =(const Polynomial&);
    /// This overloaded unary negation operator multiplies every element of the Polynomial::terms vector by -1 and then calls the simplify() method.
    Polynomial& operator -(const Polynomial&);
    /// The standard copy constructor - it calls the clear() method and then copies over all properties from the source instance to this one.
    Polynomial(const Polynomial&);
    /// This method evaluates the polynomial at the method's argument \f$\alpha\f$, returning the value of \f$p(\alpha)\f$.
    kind evaluate(kind);
    /// This method returns the value of the degree.
    unsigned int get_degree() const;
    /// This method returns true if this is the zero polynomial and false otherwise.
    bool is_null() const;
    /// This method returns the value of the coefficient specified by the method's unique argument.
    kind get_value(unsigned int) const;
    /// This method sets the coefficient specified by the second argument to the value specified by the first argument.
    void set_value(kind,unsigned int);
    /// This method clears the vector Polynomial::terms and restores all the scalar properties to their default value.
    void clear();
    /// This method writes the instance properties to a binary disk file and returns the number of bytes written to the file.
    int serialize(std::ofstream&) const;
    /// This method calls the clear() method on the instance and then reads the properties from a binary disk file and returns the number of bytes read.
    int deserialize(std::ifstream&);
    /// This method computes and returns the formal derivative of the polynomial, i.e. \f$p'(x) = d a_d x^{d-1} + \dots 2 a_2 x + a_1\f$.
    Polynomial<kind> derivative() const;
    /// This overloading of the ostream operator writes the polynomial to the screen in a "pretty print" format.
    friend std::ostream& operator << <>(std::ostream&,const Polynomial<kind>&);
    /// This overloaded operator tests the two polynomial for equality, first checking if they have the same degree and then testing the two coefficient vectors element by element.
    friend bool operator == <>(const Polynomial<kind>&,const Polynomial<kind>&);
    /// This overloaded addition operator adds together the two polynomial arguments according to the standard mathematical convention and then calls the simplify() method on the resulting output.
    friend Polynomial<kind> operator +<>(const Polynomial<kind>&,const Polynomial<kind>&);
    /// This overloaded multiplication operator multiplies a polynomial (the second argument) by a scalar (the first argument) and then calls the simplify() method on the resulting output.
    friend Polynomial<kind> operator *<>(kind,const Polynomial<kind>&);
    /// This overloaded multiplication operator multiplies two polynomials together according to the standard mathematical convention and then calls the simplify() method on the resulting output.
    friend Polynomial<kind> operator *<>(const Polynomial<kind>&,const Polynomial<kind>&);
  };

  template<class kind>
  inline unsigned int Polynomial<kind>::get_degree() const
  {
    return degree;
  }

  template<>
  /// This method is an instantiation of is_null() for the case of a polynomial over the rationals, needed to handle this case where the coefficients can be tested if they are exactly zero.
  inline bool Polynomial<Rational>::is_null() const
  {
    if (terms.empty()) return true;
    unsigned int i;
    for(i=0; i<=degree; ++i) {
      if (!terms[i].is_null()) return false;
    }
    return true;
  }

  template<class kind> 
  inline bool Polynomial<kind>::is_null() const
  {
    if (terms.empty()) return true;
    unsigned int i;
    for(i=0; i<=degree; ++i) {
      if (std::abs(terms[i]) > std::numeric_limits<double>::epsilon()) return false;
    }
    return true;
  }

  template<class kind>
  std::ostream& operator <<(std::ostream& s,const Polynomial<kind>& source)
  {
    unsigned int i;

    if (source.terms[source.degree] > 0) {
      if (source.monic) {
        if (source.degree == 1) {
          s << "x ";
        }
        else {
          s << "x^" << source.degree << " ";
        }
      }
      else {
        if (source.degree == 1) {
          s << source.terms[1] << "*x ";
        }
        else {
          s << source.terms[source.degree] << "*x^" << source.degree << " ";
        }
      }
    }
    else {
      if (source.terms[source.degree] == -1) {
        if (source.degree == 1) {
          s << "-x ";
        }
        else {
          s << "-x^" << source.degree << " ";
        }
      }
      else {
        if (source.degree == 1) {
          s << "-" << -source.terms[1] << "*x ";
        }
        else {
          s << "-" << -source.terms[source.degree] << "*x^" << source.degree << " "; 	
        }
      }
    }
    for(i=source.degree-1; i>=1; --i) {
      if (source.terms[i] == 0) continue;
      if (i > 1) {
        if (source.terms[i] == 1) {
          s << "x^" << i << " ";
        }
        else if (source.terms[i] == -1) {
          s << "-x^" << i << " ";
        }
        else {
          if (source.terms[i] > 0) {
            s << "+ " << source.terms[i] << "*x^" << i << " ";
          }
          else {
            s << "- " << -source.terms[i] << "*x^" << i << " ";
          }	    	
        }
      }	
      else if (i == 1) {
        if (source.terms[i] == 1) {
          s << "+ x ";
        }
        else if (source.terms[i] == -1) {
          s << "- x ";
        }
        else {
          if (source.terms[i] > 0) {
            s << "+ " << source.terms[i] << "*x ";
          }
          else {
            s << "- " << -source.terms[i] << "*x ";
          }	    	
        }  	
      }   	
    }
    if (!source.homogeneous) {
      if (source.terms[0] > 0) {
        s << "+ " << source.terms[0];
      }
      else {
        s << "- " << -source.terms[0]; 		
      }
    } 	
    return s;
  }

  template<>
  bool operator ==(const Polynomial<Rational>& p1,const Polynomial<Rational>& p2)
  {
    if (p1.degree != p2.degree) return false;
    unsigned int i;
    Rational q;
    for(i=0; i<p1.degree; ++i) {
      q = p1.terms[i] - p2.terms[i];
      if (!q.is_null()) return false;
    }
    return true;
  }

  template<class kind>
  bool operator ==(const Polynomial<kind>& p1,const Polynomial<kind>& p2)
  {
    if (p1.degree != p2.degree) return false;
    unsigned int i;
    for(i=0; i<p1.degree; ++i) {
      if (std::abs(p1.terms[i] - p2.terms[i]) > std::numeric_limits<double>::epsilon()) return false;
    }
    return true;
  }

  template<class kind>
  Polynomial<kind> operator -(const Polynomial<kind>& p1,const Polynomial<kind>& p2)
  {
    unsigned int i,mu = std::min(p1.degree,p2.degree);
    std::vector<kind> new_terms;

    for(i=0; i<=mu; ++i) {
      new_terms.push_back(p1.terms[i] - p2.terms[i]);
    }
    if (p1.degree > p2.degree) {
      for(i=1+mu; i<=p1.degree; ++i) {
        new_terms.push_back(p1.terms[i]);
      }
    }
    else if (p2.degree > p1.degree) {
      for(i=1+mu; i<=p2.degree; ++i) {
        new_terms.push_back(-p2.terms[i]);
      }
    }
    Polynomial<kind> output(new_terms);
    output.simplify();
    return output;
  }

  template<class kind>
  Polynomial<kind> operator +(const Polynomial<kind>& p1,const Polynomial<kind>& p2)
  {
    unsigned int i,mu = std::min(p1.degree,p2.degree);
    std::vector<kind> new_terms;

    for(i=0; i<=mu; ++i) {
      new_terms.push_back(p1.terms[i] + p2.terms[i]);
    }
    if (p1.degree > p2.degree) {
      for(i=1+mu; i<=p1.degree; ++i) {
        new_terms.push_back(p1.terms[i]);
      }
    }
    else if (p2.degree > p1.degree) {
      for(i=1+mu; i<=p2.degree; ++i) {
        new_terms.push_back(p2.terms[i]);
      }
    }
    Polynomial<kind> output(new_terms);
    output.simplify();
    return output;
  }

  template<class kind>
  Polynomial<kind> operator *(kind alpha,const Polynomial<kind>& p)
  {
    unsigned int i;
    std::vector<kind> new_terms;

    for(i=0; i<=p.degree; ++i) {
      new_terms.push_back(alpha*p.terms[i]);
    }
    Polynomial<kind> output(new_terms);
    output.simplify();
    return output;
  }

  template<class kind>
  Polynomial<kind> operator *(const Polynomial<kind>& p1,const Polynomial<kind>& p2)  
  {
    unsigned int i,j,k,mdegree = p1.degree + p2.degree;
    kind sum;
    std::vector<kind> new_terms;

    new_terms.push_back(p1.terms[0]*p2.terms[0]);
    for(i=0; i<mdegree-1; ++i) {
      // Need all the two factor additive partitions of "i"
      sum = 0;
      for(j=0; j<=p1.degree; ++j) {
        for(k=0; k<=p2.degree; ++k) {
          if ((j+k) == i) sum = sum + p1.terms[j]*p2.terms[k];
        }
      }
      new_terms.push_back(sum);
    }    
    new_terms.push_back(p1.terms[p1.degree]*p2.terms[p2.degree]);        
    Polynomial<kind> output(new_terms);
    output.simplify();
    return output;
  } 
} 
#endif





