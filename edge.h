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

#ifndef _edgeh
#define _edgeh

namespace SYNARMOSMA {
  class Edge {
   protected:
    double length;
    bool cyclic;
    double flow;
    double capacity;
    RELATION direction;
    std::set<int> vertices;

    void clear();
    void serialize(std::ofstream&) const;
    void deserialize(std::ifstream&);
   public:
    Edge();
    Edge(int,int);
    Edge(int,int,RELATION);
    Edge(const Edge&);
    ~Edge();
    Edge& operator =(const Edge&);
    inline void get_vertices(int*) const;
    inline void set_vertices(int u,int v);
    inline RELATION get_direction() const {return direction;};
    inline void set_direction(RELATION rho) {direction = rho;};
    friend class Graph;
    friend class Directed_Graph;
    friend class Logic_Graph;
  };

  void Edge::set_vertices(int u,int v) 
  {
#ifdef DEBUG
    assert(u != v);
#endif
    vertices.clear();
    vertices.insert(u);
    vertices.insert(v);
  }

  void Edge::get_vertices(int* vx) const
  {
    int i = 0;
    std::set<int>::const_iterator it;
    for(it=vertices.begin(); it!=vertices.end(); ++it) {
      vx[i] = *it; ++i;
    }
  }
}
#endif

