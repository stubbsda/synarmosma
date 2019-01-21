#include "random.h"
#include "rational.h"

#ifndef _varietyh
#define _varietyh

namespace SYNARMOSMA {
  template<class kind>
  class Variety;

  template<class kind>
  std::ostream& operator <<(std::ostream&,const Variety<kind>&);

  template<class kind>
  /// A class representing an algebraic variety, i.e. a finite system of algebraic equations in a set of unknowns, using the Monomial class. 
  class Variety {
   protected:
    /// This non-negative integer property stores the number of equations in 
    /// the variety.
    unsigned int nequation = 0;
    /// This non-negative integer property stores the number of variables in 
    /// the variety and which are enumerated successively from zero, so \f$x_0, x_1, \dots, x_{n-1}\f$.
    unsigned int nvariable = 0;
    /// This non-negative integer property is the characteristic of the domain 
    /// over which this variety is defined: zero for the integers and rationals, 
    /// a prime number p for the Galois field GF(p).
    unsigned int characteristic = 0;
    /// This Boolean property is true if all of the equations in the variety are 
    /// linear and false otherwise.
    bool linear = false;
    /// This Boolean property is true if all the equations have a remainder term of 
    /// zero and false otherwise.
    bool homogeneous = false;
    /// This Boolean property is true if every term in a given equation has the same 
    /// total degree and false otherwise.
    bool projective = false;
    std::vector<Monomial<kind> >* equations;
    std::vector<kind> remainder;
    /// This property is a vector of integer sets, each element of which contains the 
    /// independent variables upon which this equation in the variety depends.
    std::vector<std::set<unsigned int> > dependencies;
    static const kind zero;      

    /// This method allocates the memory for the array Variety::equations.
    void allocate();
    void initialize();
    void compute_properties();
    void find_partial(std::vector<unsigned int>&,int,const std::vector<unsigned int>*) const;
    int write_equations(std::ofstream&) const;
    int read_equations(std::ifstream&);
    kind generate_coefficient(int) const;
   public:
    Variety();
    Variety(int);
    Variety(int,int);
    Variety(const Variety&);
    Variety& operator =(const Variety&);
    ~Variety();
    int serialize(std::ofstream&) const;
    int deserialize(std::ifstream&);
    void elaborate();
    bool add_term(int,kind,const std::vector<unsigned int>&);
    bool add_term(int,const Monomial<kind>&);
    inline void set_remainder(int,kind);
    void make_projective();
    void clear();
    int compute_dependencies(std::vector<unsigned int>&) const;
    friend std::ostream& operator << <>(std::ostream&,const Variety<kind>&);
  };

  template<class kind>
  void Variety<kind>::set_remainder(int n,kind r)
  {
    if (n < 0 || n >= (signed) nequation) throw std::invalid_argument("Illegal equation number in Variety::set_remainder!");

    remainder[n] = r;
  }

  template<class kind>
  std::ostream& operator <<(std::ostream& s,const Variety<kind>& source)
  {
    int i;
    unsigned int j,k;
    Monomial<kind> term;

    for(i=0; i<source.nequation; ++i) {
      for(j=0; j<source.equations[i].size(); ++j) {
        term = source.equations[i][j];
        if (term.coefficient != kind(1)) s << term.coefficient << "*";
        for(k=0; k<term.exponents.size()-1; ++k) {
          s << "x(" << term.exponents[k].first << ")";
          if (term.exponents[k].second > 1) s << "^" << term.exponents[k].second;
          s << "*";
        }
        s << "x(" << term.exponents[term.exponents.size()-1].first << ")^" << term.exponents[term.exponents.size()-1].second;
        if (j < source.equations[i].size()-1) s << " + ";
      }
      if (source.remainder[i] > kind(0)) s << " + " << source.remainder[i];
      s << " = 0" << std::endl;
    }
    return s;
  }
}
#endif
