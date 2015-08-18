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

#ifndef _schemah
#define _schemah

namespace SYNARMOSMA {
  class Schema {
   protected:
    int nvertex;
    std::vector<std::set<int> > neighbours;

   public:
    Schema();
    Schema(int);
    virtual ~Schema();
    bool connected() const;
    bool connected(int,int) const;
    virtual void clear();
    int add_vertex();
    bool add_edge(int,int);
    bool drop_edge(int,int);
    bool positive_valence() const;
    int spanning_tree(std::vector<int>&) const;
    int component_analysis(std::vector<int>&) const;
    bool consistent() const;
    void components(std::vector<int>&,std::vector<int>&) const;
    inline int get_order() const {return nvertex;};
  };
}
#endif
