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

#ifndef _wordh
#define _wordh

namespace SYNARMOSMA {
  class Word {
   private:
    std::vector<std::pair<unsigned int,int> > content;
    // The number of letters in the alphabet for words in a group presentation
    unsigned int NL;

    void initialize(unsigned int);
    void initialize(const std::vector<unsigned int>&,const std::vector<int>&);
    void clear();
   public:
    Word();
    Word(unsigned int);
    Word(unsigned int,unsigned int);
    Word(unsigned int,unsigned int,int);
    Word(const Word&);
    Word& operator =(const Word&);
    ~Word();
    Word operator !() const;
    Word invert() const;
    Word mutate() const;
    Word normalize() const;
    Word swap(unsigned int,unsigned int,bool) const;
    Word reduce(int,const std::set<unsigned int>&,const unsigned int*) const;
    unsigned int size() const;
    void permute(unsigned int,Word&) const;
    bool trivial() const;
    bool alias() const;
    bool legal() const;
    bool empty() const;
    void serialize(std::ofstream&) const;
    void deserialize(std::ifstream&);
    void initialize(unsigned int,unsigned int,int);
    friend bool operator ==(const Word&,const Word&);
    friend bool operator !=(const Word&,const Word&);
    friend Word operator *(const Word&,const Word&);
    friend std::ostream& operator <<(std::ostream&,const Word&);
    friend class Group;
    friend class Homotopy;
  };
}
#endif
