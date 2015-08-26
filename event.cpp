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
  incept = -1;
  colour = -1;
  energy = 0.0;
  topological_dimension = 0;
  past.clear();
  future.clear();
  entourage.clear();
  neighbours.clear();
  proper_time.clear();
  observation.clear();
}

void Event::serialize(std::ofstream& s) const
{
  int i,j,n,m;
  char q;
  std::set<int>::const_iterator it;

  s.write((char*)(&incept),sizeof(int));
  s.write((char*)(&colour),sizeof(int));
  s.write((char*)(&topological_dimension),sizeof(int));
  s.write((char*)(&energy),sizeof(double));
  proper_time.serialize(s);
  observation.serialize(s);
  n = (signed) past.size();
  s.write((char*)(&n),sizeof(int));
  for(it=past.begin(); it!=past.end(); ++it) {
    i = *it;
    s.write((char*)(&i),sizeof(int));
  }
  n = (signed) future.size();
  s.write((char*)(&n),sizeof(int));
  for(it=future.begin(); it!=future.end(); ++it) {
    i = *it;
    s.write((char*)(&i),sizeof(int));
  }
  n = (signed) neighbours.size();
  s.write((char*)(&n),sizeof(int));
  for(it=neighbours.begin(); it!=neighbours.end(); ++it) {
    i = *it;
    s.write((char*)(&i),sizeof(int));
  }
  n = (signed) entourage.size();
  s.write((char*)(&n),sizeof(int));
  for(i=0; i<n; ++i) {
    m = (signed) entourage[i].size();
    s.write((char*)(&m),sizeof(int));
    for(j=0; j<m; ++j) {
      q = entourage[i][j];
      s.write((char*)(&q),sizeof(char));
    }
  }
}

void Event::deserialize(std::ifstream& s)
{
  int i,j,n,m;
  char q;
  std::string estring;

  clear();

  s.read((char*)(&incept),sizeof(int));
  s.read((char*)(&colour),sizeof(int));
  s.read((char*)(&topological_dimension),sizeof(int));
  s.read((char*)(&energy),sizeof(double));
  proper_time.deserialize(s);
  observation.deserialize(s);
  s.read((char*)(&n),sizeof(int));
  for(i=0; i<n; ++i) {
    s.read((char*)(&j),sizeof(int));
    past.insert(j);
  }
  s.read((char*)(&n),sizeof(int));
  for(i=0; i<n; ++i) {
    s.read((char*)(&j),sizeof(int));
    future.insert(j);
  }
  s.read((char*)(&n),sizeof(int));
  for(i=0; i<n; ++i) {
    s.read((char*)(&j),sizeof(int));
    neighbours.insert(j);
  }
  s.read((char*)(&n),sizeof(int));
  for(i=0; i<n; ++i) {
    s.read((char*)(&m),sizeof(int));
    for(j=0; j<m; ++j) {
      s.read((char*)(&q),sizeof(char));
      estring += q;
    }
    entourage.push_back(estring);
    estring.clear();
  }
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









