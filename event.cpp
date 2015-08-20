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
  initialize();
}

Event::Event(const Event& source)
{
  proper_time = source.proper_time;
  energy = source.energy;
  neighbours = source.neighbours;
  entourage = source.entourage;
  observation = source.observation;
  incept = source.incept;
  topological_dimension = source.topological_dimension;
  colour = source.colour;
  past = source.past;
  future = source.future;
}

Event& Event::operator =(const Event& source)
{
  if (this == &source)  return *this;
  proper_time = source.proper_time;
  energy = source.energy;
  neighbours = source.neighbours;
  entourage = source.entourage;
  observation = source.observation;
  incept = source.incept;
  topological_dimension = source.topological_dimension;
  colour = source.colour;
  past = source.past;
  future = source.future;
  return *this;
}

Event::~Event()
{

}

void Event::clear()
{
  initialize();
}

void Event::initialize()
{
  topological_dimension = 1;
  incept = 0;
  colour = 2;
  proper_time.initialize(RND.nrandom(1.5));
}

double Event::norm() const
{
  double output = double(topological_dimension)*proper_time.norm();
  return output;
}









