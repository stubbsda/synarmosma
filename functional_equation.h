#include "variety.h"
#include "integer_polynomial.h"

#ifndef _functionaleqnh
#define _functionaleqnh

namespace SYNARMOSMA {
  template<class kind>
  class Functional_Equation;

  extern template const int Integer_Polynomial<int>::zero;
  extern template const NTL::ZZ Integer_Polynomial<NTL::ZZ>::zero;

  template<class kind>
  std::ostream& operator <<(std::ostream& s,const Functional_Equation<kind>&);

  template<class kind>
  /// A template class representing a polynomial functional equation over a base type, which is the field for the equation's polynomial coefficients and arguments.
  class Functional_Equation {
   protected:
    /// This property contains the individual terms of the functional equation, in the form of a triple, which 
    /// represent the mathematical expression \f$\sum_{i=1}^N q_i(x) F^i(p_i(x))\f$, where \f$q_i, p_i\in K[x]\f$ 
    /// for all \f$1\le i\le N\f$, with \f$K\f$ the base field of this template class. We use the Polynomial class 
    /// to store the \f$p_i(x)\f$ and \f$q_i(x)\f$.
    std::vector<std::tuple<Integer_Polynomial<kind>,Integer_Polynomial<kind>,unsigned int> > terms;
    /// This property stores the inhomogeneous term \f$q_0(x)\in K[x]\f$ in the functional equation \f$\sum_{i=1}^N 
    /// q_i(x) F^i(p_i(x)) + q_0(x) = 0\f$, assuming it exists. 
    Integer_Polynomial<kind> remainder;
    /// This Boolean property is true if the functional equation has the form \f$q_1(x) F(p_1(x)) + q_0(x) = 0\f$, 
    /// where \f$p_1\f$ and the \f$q_i\f$ are members of \f$K[x]\f$ with \f$K\f$ the base field of this template class. 
    bool linear = false;
    /// This Boolean property is true if the functional equation has the form \f$\sum_{i=1}^N q_i(x) F^i(p_i(x)) = 0\f$,
    /// where \f$q_i, p_i\in K[x]\f$ for all \f$1\le i\le N\f$, with \f$K\f$ the base field of this template class.
    bool homogeneous = false;
   
    /// This method initializes the functional equation, with the argument representing the degree \f$N\f$ of the equation \f$\sum_{i=1}^N q_i(x) F^i(p_i(x)) + q_0(x)=0\f$. The \f$p_i, q_i\in K[x]\f$ are chosen at random and have a degree that is the minimum among 5 and \f$1+2N\f$. 
    void initialize(unsigned int);
    /// This method parses three vectors of strings to build the contents of the Functional_Equation::terms and Functional_Equation::remainder properties; the three arguments contain the coefficient polynomial \f$q_i(x)\f$, the argument polynomial \f$p_i(x)\f$ and the exponent \f$i\f$.
    void analyze_file(std::vector<std::string>&,std::vector<std::string>&,std::vector<std::string>&);
    /// This method checks that the class instance is consistent with the mathematical model of a polynomial functional equation, i.e. the degree of each term in the equation is unique, and returns true if this is so.
    bool consistent() const;
    /// This method eliminates terms whose coefficient polynomial \f$q_i(x)\f$ is identically zero and correctly sets the value of the Functional_Equation::linear and Functional_Equation::homogeneous properties; it returns true if any changes have been made to the instance.
    bool simplify();
   public:
    /// The default constructor which does nothing.
    Functional_Equation();
    /// This constructor creates a random functional equation using the initialize() method, one whose degree \f$N\ge 1\f$ is given by the argument.
    Functional_Equation(unsigned int);
    /// This constructor builds a functional equation based on the contents of a text file whose name is the argument, using the analyze_file() method.
    Functional_Equation(const std::string&);
    /// The standard copy constructor which copies over the properties from the source instance. 
    Functional_Equation(const Functional_Equation&);
    /// The standard overloaded assignment operator which copies over the properties from the source instance. 
    Functional_Equation& operator =(const Functional_Equation&);
    /// The standard destructor which in this case does nothing.
    ~Functional_Equation();
    /// This method restores all of the properties of this instance of the class to their default value.
    void clear();
    /// This method writes the properties of this instance of the class to a binary disk file and returns the number of bytes written to the file.
    int serialize(std::ofstream&) const;
    /// This method calls the clear() method on the instance and then reads the properties from a binary disk file and returns the number of bytes read.
    int deserialize(std::ifstream&);
    /// This method computes the reduction of the functional equation over the Galois field \f$GF(p^k)\f$ for prime \f$p\f$ and \f$k=1\f$. As this is a field with \f$p < \infty\f$ elements the functional equation is closed and can be converted to a variety over \f$GF(p)\f$ with \f$p\f$ equations in \f$p\f$ unknowns, {\f$F(0),F(1),\dots,F(p-1)\f$}.  
    Variety<unsigned int> reduce(unsigned int);
    /// An overloaded ostream operator to do a "pretty print" of the functional equation, using the similarly overloaded operator for the Polynomial class.  
    friend std::ostream& operator << <>(std::ostream& s,const Functional_Equation<kind>&);
  };

  template<class kind>
  std::ostream& operator <<(std::ostream& s,const Functional_Equation<kind>& source)
  {
    unsigned int i;
    std::tuple<Integer_Polynomial<kind>,Integer_Polynomial<kind>,unsigned int> trio;
    for(i=0; i<source.terms.size(); ++i) {
      trio = source.terms[i];
      s << "(" << std::get<0>(trio) << ")*F(" << std::get<1>(trio) << ")^" << std::get<2>(trio) << " +" << std::endl; 
    }
    s << source.remainder << " = 0";
    return s;
  }
}
#endif

