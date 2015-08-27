/*
  Copyright 2014 Daniel Stubbs

  This file is part of Synarmosma.

  Synarmosma is free software: you can redistribute it and/or modify 
  it under the terms of the GNU General Public License as published by 
  the Free Software Foundation, either version 3 of the License, or 
  (at your option) any later version.

  Synarmosma is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with Synarmosma.  If not, see <http://www.gnu.org/licenses/>.
*/

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
    unsigned int nequation;
    unsigned int nvariable;
    unsigned int characteristic;
    bool linear;
    bool homogeneous;
    bool projective;
    // Each element of this list contains the independent variables upon 
    // which this equation in the (algebraic) variety depends
    std::vector<std::set<unsigned int> > dependencies;
  
    void allocate();
    void initialize();
    void normalize(unsigned int);
    void property_check();
    void find_partial(std::vector<unsigned int>&,unsigned int,const std::vector<unsigned int>*) const;
    int compute_zeros() const;
   public:
    Variety();
    Variety(unsigned int);
    Variety(unsigned int,unsigned int);
    Variety(const Variety&);
    Variety& operator =(const Variety&);
    ~Variety();
    void serialize(std::ofstream&) const;
    void deserialize(std::ifstream&);
    void elaborate();
    void add_term(unsigned int,kind,const std::vector<unsigned int>&);
    void add_term(unsigned int,const Monomial<kind>&);
    void set_remainder_value(unsigned int,kind);
    void make_projective();
    void clear();
    void zeta_function(unsigned int,std::vector<unsigned int>&) const;
    int compute_dependencies(std::vector<unsigned int>&) const;
    friend std::ostream& operator << <>(std::ostream&,const Variety<kind>&);
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
