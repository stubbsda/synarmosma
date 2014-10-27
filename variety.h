#include "rational.h"
#include <NTL/ZZ_pEX.h>
#include <NTL/ZZ_pXFactoring.h>

#ifndef _varietyh
#define _varietyh

namespace SYNARMOSMA {
  template<class kind>
  class Variety;

  template<class kind>
  std::ostream& operator <<(std::ostream&,const Variety<kind>&);

  template<class kind>
  class Variety {
   private:
    std::vector<Monomial<kind> >* equations;
    std::vector<kind> remainder;
    int nequation;
    int nvariable;
    int characteristic;
    bool linear;
    bool homogeneous;
    bool projective;
    // Each element of this list contains the independent variables upon 
    // which this equation in the (algebraic) variety depends
    std::vector<std::set<int> > dependencies;
  
    void allocate();
    void initialize();
    void normalize(int);
    void find_partial(bool*,int,const std::vector<int>*) const;
    int compute_zeros();
   public:
    Variety();
    Variety(int);
    Variety(int,int);
    Variety(const Variety&);
    Variety& operator = (const Variety&);
    ~Variety();
    void elaborate();
    void add_term(int,kind,const int*);
    void add_term(int,const Monomial<kind>&);
    void set_value(int,kind);
    void make_projective();
    void clear();
    void zeta_function(int,int*);
    int compute_dependencies(int*) const;
  };

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
