#include "rational.h"
#include <NTL/ZZ_pEX.h>
#include <NTL/ZZ_pXFactoring.h>

#ifndef _varietyh
#define _varietyh

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
  void write2screen() const;
  void elaborate();
  void add_term(int,kind,const int*);
  void add_term(int,const Monomial<kind>&);
  void set_value(int,kind);
  void make_projective();
  void clear();
  void zeta_function(int,int*);
  int compute_dependencies(int*) const;
};
#endif
