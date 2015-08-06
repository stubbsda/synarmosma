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

#include "event.h"
#include "nexus.h"

#ifndef __eventspaceh
#define __eventspaceh

namespace SYNARMOSMA {
  class Eventspace {
   private:
    Event* events;  
    int nevent;
    network* causal_network;
    double L;

    void initiate_events();
    double perceptual_divergence(const double*,double,const double*,const double*) const; 
    double compute_distance(int,int) const;
    void perceptual_consistency();
   public:
    Eventspace();
    Eventspace(int);
    Eventspace(const Eventspace&);
    ~Eventspace();
    Eventspace& operator=(const Eventspace&);
    void compute_distances(std::vector<double>&) const;
    void write(const std::vector<double>&,double,const char*) const;
    void build_tessellation(Nexus*) const;
  };
}
#endif 
