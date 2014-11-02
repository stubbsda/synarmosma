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

#include "proposition.h"
#include "multitime.h"

#ifndef _eventh
#define _eventh

namespace SYNARMOSMA {
  class Event {
   private:
    int timestep;
    int colour;
    int local_dimension;
    double energy;
    Proposition observation;
    std::set<int> future,past;  
    std::set<int> neighbours;
    std::vector<std::string> entourage;
    std::vector<double> space;  
    Multitime proper_time;

    double norm() const;
    void allocate();
    void initialize();
   public:
    Event();
    Event(const Event&);
    Event(const char*);
    ~Event();  
    Event& operator=(const Event&);
    friend class Eventspace;
  };

// A class which implements a Leibnizian concept of spacetime 
/*
class Event_Graph: public Graph {
 private:    
  std::vector<unsigned int> events;     

  void perceptual_consistency();
  kind perceptual_divergence(double*,double,double*,double**) const; 
  void allocate();
  void initialize();
  
 public:
  Event_Graph();
  Event_Graph(const char*);
  Event_Graph(const Event_Graph&);
  Event_Graph & operator=(const Event_Graph&);
  ~Event_Graph();
  void initiate_events();
};
*/
}
#endif 

