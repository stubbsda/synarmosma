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
    std::vector<std::pair<bool,std::set<unsigned int> > > vertices;

   public:
    Schema();
    Schema(unsigned int);
    Schema(const Schema&);
    Schema& operator =(const Schema&);
    virtual ~Schema();
    virtual void serialize(std::ofstream&) const;
    virtual void deserialize(std::ifstream&);
    bool connected() const;
    bool connected(unsigned int,unsigned int) const;
    virtual void clear();
    unsigned int add_vertex();
    virtual bool drop_vertex(unsigned int);
    virtual bool add_edge(unsigned int,unsigned int);
    virtual bool drop_edge(unsigned int,unsigned int);
    bool positive_valence() const;
    unsigned int spanning_tree(std::vector<unsigned int>&) const;
    unsigned int component_analysis(std::vector<unsigned int>&) const;
    bool consistent() const;
    void components(std::vector<unsigned int>&,std::vector<unsigned int>&) const;
    inline unsigned int get_order() const;
  };

  inline unsigned int Schema::get_order() const 
  {
    unsigned int order = 0;
    std::vector<std::pair<bool,std::set<unsigned int> > >::const_iterator vit;

    for(vit=vertices.begin(); vit!=vertices.end(); ++vit) {
      if (vit->first) order++;
    }
    return order;
  }
}
#endif
