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

#include "cell.h"
#include "graph.h"

#ifndef _nexush
#define _nexush

namespace SYNARMOSMA {
  class Nexus : public Schema {
   private:
    int dimension;
    std::vector<Cell>* elements;
    hash_map* index_table;

   public:
    Nexus();
    Nexus(int);
    virtual ~Nexus();
    bool orientable() const;
    bool pseudomanifold(bool*) const;
    void assemble();
    void surface_construction(int);
    virtual void clear();
    void initialize(int);
    void initialize(int,int);
    void paste(const std::set<int>&);
    void regularization();
    int size() const;
    void ascend(int,int,std::vector<Cell>&) const;
    void star(const std::set<std::set<int> >&,std::vector<Cell>*) const;
    void link(const std::set<std::set<int> >&,std::vector<Cell>*) const;
    void closure(const std::set<std::set<int> >&,Nexus*,int*) const;
    void compute_entourages();
    void compute_neighbours();
    inline int get_dimension() const {return dimension;};
    friend class Homology;
    friend class Homotopy;
  };
}
#endif


