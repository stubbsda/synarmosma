#include "variety.h"
#include "polynomial.h"

#ifndef _functionaleqnh
#define _functionaleqnh

namespace SYNARMOSMA {
  /// A class representing a polynomial functional equation, where the base field for the equation's polynomial coefficients and arguments is the rationals.
  class Functional_Equation {
   protected:
    /// This property contains the individual terms of the functional equation, in the form of a triple, which
    /// represent the mathematical expression \f$\sum_{i=1}^N \alpha_i(x) [F(\beta_i(x))]^i\f$, where the
    /// \f$\alpha_i, \beta_i\in \mathbf{Q}[x]\f$ for \f$1\le i\le N\f$. We use the Polynomial class over 
    /// the base type Rational to store the \f$\alpha_i(x)\f$ and \f$\beta_i(x)\f$.
    std::vector<std::tuple<Polynomial<Rational>,Polynomial<Rational>,unsigned int> > terms;
    /// This property stores the inhomogeneous term \f$\gamma(x)\in \mathbf{Q}[x]\f$ in the functional equation \f$\sum_{i=1}^N
    /// \alpha_i(x) [F(\beta_i(x))]^i + \gamma(x) = 0\f$, assuming it exists.
    Polynomial<Rational> constant;
    /// This Boolean property is true if the functional equation has the form \f$\alpha_1(x) F(\beta_1(x)) + \gamma(x) = 0\f$,
    /// where \f$\alpha_1\f$, \f$\beta_1\f$ and \f$\gamma\f$ are all members of \f$\mathbf{Q}[x]\f$.
    bool linear = false;
    /// This Boolean property is true if the functional equation has the form \f$\sum_{i=1}^N \alpha_i(x) [F(\beta_i(x))]^i = 0\f$,
    /// where \f$\alpha_i, \beta_i\in \mathbf{Q}[x]\f$ for all \f$1\le i\le N\f$, with the constant identically equal to zero.
    bool homogeneous = false;

    /// This method takes as its first argument a string of the form "(-1,0,2/3,7)", representing the rational polynomial \f$7x^3 + (2/3)x^2 -1\f$ in this case, and parses the string to form a vector of elements of the Rational class, written to the method's second argument, which can be used to initialize an instance of the Polynomial class.
    void parse_polynomial(const std::string&,std::vector<Rational>&) const;
    /// This method parses three vectors of strings to build the contents of the Functional_Equation::terms and Functional_Equation::constant properties; the three arguments contain the coefficient polynomial \f$\alpha_i(x)\f$, the argument polynomial \f$\beta_i(x)\f$ and the exponent \f$i\ge 1\f$ as well as the constant term \f$\gamma(x)\f$ when \f$i=0\f$. Each of these strings encoding a rational polynomial are then analyzed using parse_polynomial().
    void parse_equation(const std::vector<std::string>&,const std::vector<std::string>&,const std::vector<std::string>&);
    /// This method checks that the class instance is consistent with the mathematical model of a polynomial functional equation, i.e. the degree of each term in the equation is unique, and returns true if this is so.
    bool consistent() const;
    /// This method eliminates terms whose coefficient polynomial \f$\alpha_i(x)\f$ is identically zero and correctly sets the value of the Functional_Equation::linear and Functional_Equation::homogeneous properties; it returns true if any changes have been made to the instance.
    bool simplify();
   public:
    /// The default constructor which does nothing.
    Functional_Equation();
    /// This constructor builds a functional equation based on the contents of a text file whose name is the argument, using the parse_equation() method.
    Functional_Equation(const std::string&);
    /// The standard copy constructor which copies over the properties from the source instance. 
    Functional_Equation(const Functional_Equation&);
    /// The standard overloaded assignment operator which copies over the properties from the source instance. 
    Functional_Equation& operator =(const Functional_Equation&);
    /// The standard destructor which in this case does nothing.
    virtual ~Functional_Equation();
    /// This method restores all of the properties of this instance of the class to their default value.
    void clear();
    /// This method writes the properties of this instance of the class to a binary disk file and returns the number of bytes written to the file.
    int serialize(std::ofstream&) const;
    /// This method calls the clear() method on the instance and then reads the properties from a binary disk file and returns the number of bytes read.
    int deserialize(std::ifstream&);
    /// This method computes the reduction of the functional equation over the Galois field \f$GF(p^k)\f$ for prime \f$p\f$ and \f$k=1\f$. As this is a field with \f$p < \infty\f$ elements the functional equation is closed and can be converted to a variety over \f$GF(p)\f$ with \f$p\f$ equations in \f$p\f$ unknowns, {\f$F(0),F(1),\dots,F(p-1)\f$}.  
    bool reduce(unsigned int,Variety<unsigned int>*) const;
    /// An overloaded ostream operator to do a "pretty print" of the functional equation, using the similarly overloaded operator for the Integer_Polynomial class.  
    friend std::ostream& operator <<(std::ostream& s,const Functional_Equation&);
  };

  template<class kind>
  std::ostream& operator <<(std::ostream& s,const Functional_Equation& source)
  {
    unsigned int i;
    std::tuple<Polynomial<Rational>,Polynomial<Rational>,unsigned int> trio;
    for(i=0; i<source.terms.size(); ++i) {
      trio = source.terms[i];
      s << "(" << std::get<0>(trio) << ")*F(" << std::get<1>(trio) << ")^" << std::get<2>(trio) << " +" << std::endl; 
    }
    s << source.constant << " = 0";
    return s;
  }
}
#endif

