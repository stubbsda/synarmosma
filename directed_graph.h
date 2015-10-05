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

#include "poset.h"
#include "graph.h"

#ifndef __dgraph
#define __dgraph

namespace SYNARMOSMA {
  class Directed_Graph : public Graph {
   private:
    // Hyphantic operators
    bool foliation_x(int,int);
    bool foliation_m(int,int);
    int fission_x(int);
    int fission_m(int);
   public:
    Directed_Graph();
    Directed_Graph(int);
    Directed_Graph(int,double);
    Directed_Graph(const Directed_Graph&);
    Directed_Graph& operator =(const Directed_Graph&);
    virtual ~Directed_Graph();
    virtual void clear();
    virtual void serialize(std::ofstream&) const;
    virtual void deserialize(std::ifstream&);
    virtual bool consistent() const;
    virtual int add_vertex();
    virtual bool drop_vertex(int);
    virtual bool drop_edge(int,int);
    virtual bool add_edge(int,int);
    bool add_edge(int,int,RELATION);
    int directedness() const;
    bool path_connected(int,int) const;
    void compute_sinks(std::set<int>&) const;
    void compute_sources(std::set<int>&) const;
  };
}
#endif
