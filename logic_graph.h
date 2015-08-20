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
#include "propositional_system.h"

#ifndef _lgraphh
#define _lgraphh

namespace SYNARMOSMA {
  class Logic_Graph: public Graph {
   private:  
    // The main property of this class
    Propositional_System* logic;
    // This array measures the number of atomic propositions used in a given 
    // graph neighbourhood (a vertex and its neighbours)
    std::vector<int> logical_breadth;

    double rationalize_topology();
    void compute_logical_breadth();

    virtual bool amputation(int);
    virtual bool fusion(int,int);
    virtual bool foliation_x(int,int);
    virtual bool foliation_m(int,int);
    virtual int fission_x(int);
    virtual int fission_m(int);
    virtual bool add_edge(int,int);  
   public:
    Logic_Graph();
    Logic_Graph(int);
    Logic_Graph(const Logic_Graph&);
    Logic_Graph& operator =(const Logic_Graph&);
    virtual ~Logic_Graph();
    void create();
  };
}
#endif
