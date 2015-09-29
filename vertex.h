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

#ifndef _eventh
#define _eventh

namespace SYNARMOSMA {
  class Vertex {
   protected:
    int incept;
    int topological_dimension;
    double energy;
    Proposition theorem;
    std::set<int> posterior,anterior;  
    std::set<int> neighbours;
    std::set<int> entourage;

   public:
    Vertex();
    Vertex(const Vertex&);
    ~Vertex();  
    Vertex& operator =(const Vertex&);
    void serialize(std::ofstream&) const;
    void deserialize(std::ifstream&);
    void clear();
  };
}
#endif 

