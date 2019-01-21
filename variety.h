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
    /// This vector stores the remainder term for each equation in the variety and 
    /// thus should have a length equal to Variety::nequation. 
    std::vector<kind> remainder;
    /// This property is a vector of integer sets, each element of which contains the 
    /// independent variables upon which this equation in the variety depends.
    std::vector<std::set<unsigned int> > dependencies;
    static const kind zero;      

    /// This method allocates the memory for the array Variety::equations.
    void allocate();
    /// This method tries to determine the correct value of various of the class properties, such as Variety::homogeneous, Variety::linear and Variety::projective. 
    void compute_properties();
    void find_partial(std::vector<unsigned int>&,int,const std::vector<unsigned int>*) const;
    int write_equations(std::ofstream&) const;
    int read_equations(std::ifstream&);
    /// This utility method generates a random number from the appropriate base field between -L and L, where L is the method's argument, while also respecting the field charactertistic. 
    kind generate_coefficient(int) const;
   public:
    /// The default constructor which does nothing.
    Variety();
    /// The principal constructor for the Variety class, which sets the Variety::nequation property to the argument and then calls the allocate() method.
    Variety(int);
    /// This constructor sets Variety::nequation equal to the first argument and the set argument is the characteristic of the base field; the constructor then calls the allocate() method.
    Variety(int,int);
    /// The standard copy constructor which copies over the properties from the source instance. 
    Variety(const Variety&);
    /// The standard overloaded assignment operator which copies over the properties from the source instance. 
    Variety& operator =(const Variety&);
    /// The destructor which, if Variety::nequation is greater than -1, frees the memory in the Variety::equations property.
    ~Variety();
    /// This method writes the properties of this instance of the class to a binary disk file and returns the number of bytes written to the file.
    int serialize(std::ofstream&) const;
    /// This method calls the clear() method on the instance and then reads the properties from a binary disk file and returns the number of bytes read.
    int deserialize(std::ifstream&);
    void elaborate();
    /// This method initializes a random algebraic variety using the generate_coefficient() method for the coefficients and remainder term, then calls compute_properties(). The method's argument is the maximum total degree of each term in each equation. 
    void random_variety(int);
    /// This method adds a term to the equation specified by the first argument; the second argument is the term's coefficient and the third is the vector of non-negative powers for each variable in the variety. 
    bool add_term(int,kind,const std::vector<unsigned int>&);
    bool add_term(int,const Monomial<kind>&);
    /// This method sets the remainder to the value in the second argument for the equation specified by the first argument. 
    inline void set_remainder(int,kind);
    void make_projective();
    /// This method frees the memory associated with the Variety::equations property (if any has been allocated), clears the Variety::remainder and Variety::dependencies vectors and sets all of the class' properties back to their default values.
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
