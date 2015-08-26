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

#include "proposition.h"
#include "binary_matrix.h"
#include "directed_graph.h"

#ifndef _psystemh
#define _psystemh

namespace SYNARMOSMA {
  class Propositional_System {
   private:
    std::vector<Proposition> theorems;
    unsigned int natom;
    unsigned int nuniverse;
    std::vector<boost::dynamic_bitset<> > truth;

    void compute_internal_logic();
    void set_default_values();
    void initialize(unsigned int);
  
   public:
    Propositional_System();
    Propositional_System(unsigned int);
    Propositional_System(unsigned int,unsigned int);
    Propositional_System(unsigned int,const char*);
    Propositional_System(const Propositional_System&);
    Propositional_System& operator =(const Propositional_System&);
    ~Propositional_System();
    void clear();
    void serialize(std::ofstream&) const;
    void deserialize(std::ifstream&);
    unsigned int bit_count(unsigned int) const;
    unsigned int consistency(unsigned int,unsigned int,const std::string&) const;
    bool implication(unsigned int,const std::vector<unsigned int>&) const;
    void compute_implication_graph(Directed_Graph*) const;
    friend class Logic_Graph;
  };
}
#endif 
