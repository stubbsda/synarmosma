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
   private:
    bool active;
    double length;
    bool cyclic;
    double flow;
    double capacity;
    std::set<int> nodes;

    void clear();
    void serialize(std::ofstream&) const;
    void deserialize(std::ifstream&);
   public:
    Edge();
    Edge(int,int);
    Edge(const Edge&);
    ~Edge();
    Edge& operator =(const Edge&);
    inline void get_nodes(int*) const;
    inline bool get_activity() const {return active;};
    inline void set_activity(bool a) {active = a;};
    friend class Graph;
    friend class Directed_Graph;
    friend class Logic_Graph;
  };

  inline void Edge::get_nodes(int* vx) const
  {
    int i = 0;
    std::set<int>::const_iterator it;
    for(it=nodes.begin(); it!=nodes.end(); ++it) {
      vx[i] = *it; ++i;
    }
  }
}
#endif
