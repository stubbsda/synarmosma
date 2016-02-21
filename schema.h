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
    Schema(const Schema&);
    Schema& operator =(const Schema&);
    virtual ~Schema();
    virtual void serialize(std::ofstream&) const;
    virtual void deserialize(std::ifstream&);
    bool connected() const;
    inline bool connected(int,int) const;
    virtual void clear();
    inline int add_vertex();
    virtual bool add_edge(int,int);
    virtual bool drop_edge(int,int);
    virtual int distance(int,int) const;
    bool positive_valence() const;
    int spanning_tree(std::vector<int>&) const;
    int component_analysis(std::vector<int>&) const;
    virtual bool consistent() const;
    void components(std::vector<int>&,std::vector<int>&) const;
    inline int get_order() const {return nvertex;};
  };

  int Schema::add_vertex()
  {
    std::set<int> empty;
    neighbours.push_back(empty);
    nvertex++;
    return nvertex-1;
  }

  bool Schema::connected(int n,int m) const
  {
    // A method to check if the vertices n and m share a direct 
    // connection
#ifdef DEBUG
    assert(n >= 0 && n < nvertex);
    assert(m >= 0 && m < nvertex);
    if (n == m) return false;
#endif

    if (neighbours[n].count(m) == 0) {
      // This edge doesn't exist...
#ifdef DEBUG
      assert(neighbours[m].count(n) == 0);
#endif
      return false;
    }
    return true;
  }
}
#endif
