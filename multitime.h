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

#ifndef _mtimeh
#define _mtimeh

namespace SYNARMOSMA {
  class Multitime {
   private:
    std::vector<std::pair<double,bool> > chronos; 
    // Not too Aristotelian, perhaps?
    static const unsigned int tdimension = 1;

    void clear();
    void allocate();
    void initialize(double);
   public:
    Multitime();
    Multitime(double);
    Multitime(const Multitime&);
    ~Multitime();
    Multitime& operator =(const Multitime&);
    void serialize(std::ofstream&) const;
    void deserialize(std::ifstream&);
    double norm() const;
    void extract(std::vector<double>&) const;
    friend Multitime operator +(const Multitime&,const Multitime&);
    friend Multitime operator -(const Multitime&,const Multitime&);
    friend Multitime operator *(const Multitime&,const Multitime&);
    friend bool operator >(const Multitime&,const Multitime&);
    friend bool operator <(const Multitime&,const Multitime&);
    friend bool operator ==(const Multitime&,const Multitime&);
    friend bool operator !=(const Multitime&,const Multitime&);
    friend bool operator >=(const Multitime&,const Multitime&);
    friend bool operator <=(const Multitime&,const Multitime&);
    friend Multitime operator *(double alpha,const Multitime&);
    friend class Event;
  };
}
#endif 
