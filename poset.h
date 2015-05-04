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
    int N;
    boost::unordered_map<std::pair<int,int>,bool> order;

   public:
    Poset();
    Poset(int);
    Poset(const Poset&);
    ~Poset();
    void clear();
    bool consistent() const;
    void add_vertex();
    bool set_order(int,int); 
    void construct_order(double);
    bool covered(int,int) const;
    int chain_number(int) const;
    int width(int,int) const;
    void write_incastrature(const std::string&) const;
    RELATION get_relation(int,int) const;
  };
}
#endif

