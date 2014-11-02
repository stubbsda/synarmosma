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

#include "global.h"

#ifndef _propositionh
#define _propositionh

namespace SYNARMOSMA {
  class Proposition {
   private:
    std::vector<int> clause;

    // The number of atomic propositions in a clause
    static const int NP = 5;
   public:
    Proposition();
    Proposition(int);
    Proposition(const std::set<int>&);
    Proposition(int,const std::set<int>&);
    Proposition(const Proposition&);
    Proposition& operator =(const Proposition&);
    Proposition& operator *(const Proposition&);
    ~Proposition();
    void initialize(int,const std::set<int>&);
    bool evaluate(const std::vector<int>&,std::vector<int>&) const;
    bool evaluate(const bool*) const;
    bool satisfiable() const;
    void atoms(std::set<int>&) const;
    inline void clear() {clause.clear();};
    void serialize(std::ofstream&) const;
    void deserialize(std::ifstream&); 
    inline static int get_clause_size() {return NP;};
    friend std::ostream& operator <<(std::ostream&,const Proposition&);
    friend class Logic_Graph;  
    friend class Propositional_System;
  };
}
#endif
