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

#include "nexus.h"

#ifndef _homologyh
#define _homologyh

namespace SYNARMOSMA {
  // How to handle the base ZZ_p ?
  enum FIELD 
  {
      INT,
      ZZ,
      GF2
  };

  enum METHOD 
  {
      GAP,
      NATIVE
  };

  class Homology {
   private:
    //std::vector<Group> sequence;
    //std::vector<std::pair<unsigned int,std::vector<unsigned int> > > sequence;
    std::vector<unsigned int> betti_number;
    std::vector<std::vector<unsigned int> > torsion;
    METHOD method;
    FIELD field;

    void compute_integral_native(const Nexus*);
    void compute_native(const Nexus*);
    void compute_gap(const Nexus*);
  
   public:
    Homology();
    Homology(FIELD,METHOD);
    Homology(const Homology&);
    Homology& operator =(const Homology&);
    ~Homology();
    std::string write() const;
    inline void set_method(METHOD m) {method = m;};
    inline void set_field(FIELD f) {field = f;};
    inline METHOD get_method() const {return method;};
    inline FIELD get_field() const {return field;};
    inline void get_betti_numbers(std::vector<unsigned int>& output) const {output = betti_number;}; 
    void initialize(FIELD,METHOD);
    void clear();
    void compute(const Nexus*);
    void serialize(std::ofstream&) const;
    void deserialize(std::ifstream&);
  };
}
#endif 
