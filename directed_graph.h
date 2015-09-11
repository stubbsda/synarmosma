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

#include "graph.h"
#include "poset.h"

#ifndef __dgraph
#define __dgraph

namespace SYNARMOSMA {
  class Directed_Graph : public Graph {
   private:
    Poset* orientation;

    // Hyphantic operators
    bool foliation_x(unsigned int);
    bool foliation_x(unsigned int,unsigned int);
    bool foliation_m(unsigned int);
    bool foliation_m(unsigned int,unsigned int);
    unsigned int fission_x(unsigned int);
    unsigned int fission_m(unsigned int);
   public:
    Directed_Graph();
    Directed_Graph(unsigned int);
    Directed_Graph(const Directed_Graph&);
    Directed_Graph& operator =(const Directed_Graph&);
    virtual ~Directed_Graph();
    virtual void clear();
    virtual void serialize(std::ofstream&) const;
    virtual void deserialize(std::ifstream&);
    bool add_edge(unsigned int,unsigned int);
    bool add_edge(unsigned int,unsigned int,RELATION);
    bool alter_orientation(unsigned int,unsigned int,RELATION);
    unsigned int two_cycles() const;
    bool path_connected(unsigned int,unsigned int) const;
    void compute_sinks(std::set<unsigned int>&) const;
    void compute_sources(std::set<unsigned int>&) const;
  };
}
#endif
