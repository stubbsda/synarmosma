#include "variety.h"
#include "polynomial.h"

#ifndef _functionaleqnh
#define _functionaleqnh

namespace SYNARMOSMA {
  template<class kind>
  class Functional_Equation;

  template<class kind>
  std::ostream& operator <<(std::ostream& s,const Functional_Equation<kind>&);

  template<class kind>
  /// A template class representing a polynomial functional equation over a base type, which is the field for the equation's polynomial coefficients and arguments.
  class Functional_Equation {
   protected:
    std::vector<std::tuple<Polynomial<kind>,Polynomial<kind>,unsigned int> > terms;
    Polynomial<kind> remainder;
    /// This Boolean property is true if the functional equation has the form \f$q_1(x) F(p_1(x)) + q_0(x) = 0\f$, 
    /// where \f$p_1\f$ and the \f$q_i\f$ are members of \f$K[x]\f$ with \f$K\f$ the base field of this template class. 
    bool linear = false;
    /// This Boolean property is true if the functional equation has the form \f$\sum_{i=1}^N q_i(x) F^i(p_i(x)) = 0\f$,
    /// where \f$q_i, p_i\in K[x]\f$ for all \f$1\le i\le N\f$, with \f$K\f$ the base field of this template class.
    bool homogeneous = false;
   
    void initialize(unsigned int);
    void analyze_file(std::vector<std::string>&,std::vector<std::string>&,std::vector<std::string>&);
   public:
    /// The default constructor which does nothing.
    Functional_Equation();
    Functional_Equation(unsigned int);
    Functional_Equation(const std::string&);
    Functional_Equation(const Functional_Equation&);
    Functional_Equation& operator =(const Functional_Equation&);
    ~Functional_Equation();
    void clear();
    /// This method writes the properties of this instance of the class to a binary disk file and returns the number of bytes written to the file.
    int serialize(std::ofstream&) const;
    /// This method calls the clear() method on the instance and then reads the properties from a binary disk file and returns the number of bytes read.
    int deserialize(std::ifstream&);
    Variety<unsigned int> reduce(unsigned int);  
    friend std::ostream& operator << <>(std::ostream& s,const Functional_Equation<kind>&);
  };

  template<class kind>
  std::ostream& operator <<(std::ostream& s,const Functional_Equation<kind>& source)
  {
    unsigned int i;
    std::tuple<Polynomial<kind>,Polynomial<kind>,unsigned int> trio;
    for(i=0; i<source.terms.size(); ++i) {
      trio = source.terms[i];
      s << "(" << std::get<0>(trio) << ")*F(" << std::get<1>(trio) << ")^" << std::get<2>(trio) << " +" << std::endl; 
    }
    s << source.remainder << " = 0";
    return s;
  }
}
#endif

