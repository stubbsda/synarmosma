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

#include "vertex.h"

using namespace SYNARMOSMA;

extern Random RND;

Vertex::Vertex()
{
  clear();
}

Vertex::Vertex(const Vertex& source)
{
  energy = source.energy;
  neighbours = source.neighbours;
  entourage = source.entourage;
  theorem = source.theorem;
  incept = source.incept;
  topological_dimension = source.topological_dimension;
  anterior = source.anterior;
  posterior = source.posterior;
}

Vertex& Vertex::operator =(const Vertex& source)
{
  if (this == &source)  return *this;
  energy = source.energy;
  neighbours = source.neighbours;
  entourage = source.entourage;
  theorem = source.theorem;
  incept = source.incept;
  topological_dimension = source.topological_dimension;
  anterior = source.anterior;
  posterior = source.posterior;
  return *this;
}

Vertex::~Vertex()
{

}

void Vertex::clear()
{
  incept = -1;
#ifdef DISCRETE
  energy = 0;
#else
  energy = 0.0;
#endif
  topological_dimension = 0;
  anterior.clear();
  posterior.clear();
  entourage.clear();
  neighbours.clear();
  theorem.clear();
}

void Vertex::serialize(std::ofstream& s) const
{
  int i,n;
  std::set<int>::const_iterator it;

  s.write((char*)(&incept),sizeof(int));
  s.write((char*)(&topological_dimension),sizeof(int));
#ifdef DISCRETE
  s.write((char*)(&energy),sizeof(UINT64));
#else
  s.write((char*)(&energy),sizeof(double));
#endif
  theorem.serialize(s);
  n = (signed) anterior.size();
  s.write((char*)(&n),sizeof(int));
  for(it=anterior.begin(); it!=anterior.end(); ++it) {
    i = *it;
    s.write((char*)(&i),sizeof(int));
  }
  n = (signed) posterior.size();
  s.write((char*)(&n),sizeof(int));
  for(it=posterior.begin(); it!=posterior.end(); ++it) {
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
  for(it=entourage.begin(); it!=entourage.end(); ++it) {
    i = *it;
    s.write((char*)(&i),sizeof(int));
  }
}

void Vertex::deserialize(std::ifstream& s)
{
  int i,j,n;

  clear();

  s.read((char*)(&incept),sizeof(int));
  s.read((char*)(&topological_dimension),sizeof(int));
#ifdef DISCRETE
  s.read((char*)(&energy),sizeof(UINT64));
#else
  s.read((char*)(&energy),sizeof(double));
#endif
  theorem.deserialize(s);
  s.read((char*)(&n),sizeof(int));
  for(i=0; i<n; ++i) {
    s.read((char*)(&j),sizeof(int));
    anterior.insert(j);
  }
  s.read((char*)(&n),sizeof(int));
  for(i=0; i<n; ++i) {
    s.read((char*)(&j),sizeof(int));
    posterior.insert(j);
  }
  s.read((char*)(&n),sizeof(int));
  for(i=0; i<n; ++i) {
    s.read((char*)(&j),sizeof(int));
    neighbours.insert(j);
  }
  s.read((char*)(&n),sizeof(int));
  for(i=0; i<n; ++i) {
    s.read((char*)(&j),sizeof(int));
    entourage.insert(j);
  }
}









