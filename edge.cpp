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

#include "edge.h"

using namespace SYNARMOSMA;

Edge::Edge()
{
  clear();
}

Edge::Edge(int x,int y) 
{
  assert(x != y);
  clear();
  nodes.insert(x);
  nodes.insert(y);
}


Edge::~Edge()
{

}

Edge::Edge(const Edge& source)
{
  active = source.active;
  length = source.length;
  flow = source.flow;
  capacity = source.capacity;
  nodes = source.nodes;
  cyclic = source.cyclic;
  direction = source.direction;
}

Edge& Edge::operator =(const Edge& source)
{
  if (this == &source) return *this;

  active = source.active;
  length = source.length;
  flow = source.flow;
  capacity = source.capacity;
  nodes = source.nodes;
  cyclic = source.cyclic;
  direction = source.direction;

  return *this;
}

void Edge::clear()
{
  active = true;
  length = 0.0;
  nodes.clear();
  capacity = 0.0;
  flow = 0.0;
  cyclic = false;
  direction = DISPARATE;
}

void Edge::serialize(std::ofstream& s) const
{
  int n;
  std::set<int>::const_iterator it;
  for(it=nodes.begin(); it!=nodes.end(); ++it) {
    n = *it;
    s.write((char*)(&n),sizeof(int));
  }
  s.write((char*)(&active),sizeof(bool));
  s.write((char*)(&cyclic),sizeof(bool));
  s.write((char*)(&length),sizeof(double));
  s.write((char*)(&flow),sizeof(double));
  s.write((char*)(&capacity),sizeof(double));
  s.write((char*)(&direction),sizeof(int));
}

void Edge::deserialize(std::ifstream& s)
{
  int n;

  clear();

  s.read((char*)(&n),sizeof(int));
  nodes.insert(n);
  s.read((char*)(&n),sizeof(int));
  nodes.insert(n);
  s.read((char*)(&active),sizeof(bool));
  s.read((char*)(&cyclic),sizeof(bool));
  s.read((char*)(&length),sizeof(double));
  s.read((char*)(&flow),sizeof(double));
  s.read((char*)(&capacity),sizeof(double));
  s.read((char*)(&direction),sizeof(int));
}

