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

#include "group.h"
#include "nexus.h"

#ifndef _homotopyh
#define _homotopyh

namespace SYNARMOSMA {
  class Homotopy {
   private:
    std::vector<Group> sequence;
    double fitness;

    void compute_fitness();
   public:
    Homotopy();
    Homotopy(int);
    Homotopy(const Homotopy&);
    ~Homotopy();
    std::string write() const;
    void mutate();
    void clear();
    void compute(const Nexus*);
    void serialize(std::ofstream&) const;
    void deserialize(std::ifstream&);
    inline double get_fitness() const {return fitness;};
    friend Homotopy operator +(const Homotopy&,const Homotopy&);
    Homotopy& operator =(const Homotopy&);  
  };
}
#endif
