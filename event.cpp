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

using namespace SYNARMOSMA;

extern Random RND;

Event::Event()
{
  allocate();
  initialize();
}

Event::Event(const char* filename)
{
  allocate();
  initialize();
}

Event::Event(const Event& source)
{
  proper_time = source.proper_time;
  space = source.space;
  energy = source.energy;
  neighbours = source.neighbours;
  entourage = source.entourage;
  observation = source.observation;
  timestep = source.timestep;
  local_dimension = source.local_dimension;
  colour = source.colour;
  past = source.past;
  future = source.future;
}

Event& Event::operator=(const Event& source)
{
  if (this != &source) {
    proper_time = source.proper_time;
    space = source.space;
    energy = source.energy;
    neighbours = source.neighbours;
    entourage = source.entourage;
    observation = source.observation;
    timestep = source.timestep;
    local_dimension = source.local_dimension;
    colour = source.colour;
    past = source.past;
    future = source.future;
  }
  return *this;
}

Event::~Event()
{

}

void Event::allocate()
{


}

void Event::initialize()
{
  // The obvious answer for the dimension of space
  local_dimension = 3;
  timestep = 0;
  colour = 2;
  proper_time.initialize(RND.nrandom(1.5));
  for(int i=0; i<local_dimension; ++i) {
    space.push_back(-5.0 + RND.nrandom(10.0));
  }
}

double Event::norm() const
{
  double output = 0.0;
  double rl = 0.0;
  for(int i=0; i<local_dimension; ++i) {
    rl += std::abs(space[i])*std::abs(space[i]);
  }
  rl = std::sqrt(rl);
  output = proper_time.norm() + rl;
  return output;
}









