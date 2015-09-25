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

#ifndef _poseth
#define _poseth

namespace SYNARMOSMA {
  enum RELATION
  {
      BEFORE,
      AFTER,
      INCOMPARABLE
  };

  class Poset {
   protected:
    unsigned int N;
    boost::unordered_map<std::pair<unsigned int,unsigned int>,bool> order;

    void compute_width(unsigned int,unsigned int,std::set<unsigned int>&) const;
    unsigned int build_chain(std::vector<unsigned int>&,unsigned int) const;
    virtual void serialize(std::ofstream&) const;
    virtual void deserialize(std::ifstream&);
   public:
    Poset();
    Poset(unsigned int);
    Poset(const Poset&);
    Poset& operator =(const Poset&);
    virtual ~Poset();
    virtual void clear();
    virtual bool consistent() const;
    inline void add_element() {N += 1;};
    bool sink(unsigned int) const;
    bool source(unsigned int) const;
    void compute_anteriority(unsigned int,std::set<unsigned int>&) const;
    void compute_posteriority(unsigned int,std::set<unsigned int>&) const;
    bool set_order(unsigned int,unsigned int); 
    bool unset_order(unsigned int,unsigned int); 
    bool invert_order(unsigned int,unsigned int);
    void construct_ordering(double);
    bool covered(unsigned int,unsigned int) const;
    unsigned int chain_number(unsigned int) const;
    void write_incastrature(const std::string&) const;
    RELATION get_order(unsigned int,unsigned int) const;
    friend class Directed_Graph;
  };
}
#endif

