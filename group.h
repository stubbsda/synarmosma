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

#include "word.h"

#ifndef _grouph
#define _grouph

namespace SYNARMOSMA {
  // A class for a combinatorial group presentation
  class Group {
   private:
    unsigned int ngenerator;
    std::vector<Word> relations;
    bool abelian;
    bool finite;
    bool solvable;
    bool free;
    bool braid;
    unsigned int cardinality;
    unsigned int rank;
    std::vector<unsigned int> torsion;

    void allocate(unsigned int);
    void compute_rank();
    Group abelianize() const;
   public:
    Group();
    Group(int);
    Group(int,const std::vector<Word>&);
    Group(int,int);
    Group(unsigned int,const std::vector<unsigned int>&);
    Group(const Group&);
    Group(const std::string&,int);
    ~Group();
    Group& operator =(const Group&);
    inline int get_rank() const {return rank;};
    unsigned int implied_generators() const;
    void initialize(unsigned int,const std::vector<Word>&);
    void initialize(unsigned int,const std::vector<unsigned int>&);
    bool equivalent(const Word&,const Word&) const;
    void reduce();
    void clear();
    void create_random();
    bool consistent() const;
    void serialize(std::ofstream&) const;
    void deserialize(std::ifstream&);
    std::string compact_form() const;
    friend std::ostream& operator <<(std::ostream&,const Group&);
    friend class Homotopy;
  };
}
#endif
