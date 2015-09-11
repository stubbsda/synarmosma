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

Edge::Edge(unsigned int x,unsigned int y) 
{
  clear();
  nodes[0] = x;
  nodes[1] = y;
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
  nodes[0] = source.nodes[0];
  nodes[1] = source.nodes[1];
  cyclic = source.cyclic;
}

Edge& Edge::operator =(const Edge& source)
{
  if (this == &source) return *this;

  active = source.active;
  length = source.length;
  flow = source.flow;
  capacity = source.capacity;
  nodes[0] = source.nodes[0];
  nodes[1] = source.nodes[1];
  cyclic = source.cyclic;

  return *this;
}

void Edge::clear()
{
  active = true;
  length = 0.0;
  nodes[0] = 0;
  nodes[1] = 0;
  capacity = 0.0;
  flow = 0.0;
  cyclic = false;
}

void Edge::serialize(std::ofstream& s) const
{
  s.write((char*)(&nodes[0]),sizeof(int));
  s.write((char*)(&nodes[1]),sizeof(int));
  s.write((char*)(&active),sizeof(bool));
  s.write((char*)(&cyclic),sizeof(bool));
  s.write((char*)(&length),sizeof(double));
  s.write((char*)(&flow),sizeof(double));
  s.write((char*)(&capacity),sizeof(double));
}

void Edge::deserialize(std::ifstream& s)
{
  clear();

  s.read((char*)(&nodes[0]),sizeof(int));
  s.read((char*)(&nodes[1]),sizeof(int));
  s.read((char*)(&active),sizeof(bool));
  s.read((char*)(&cyclic),sizeof(bool));
  s.read((char*)(&length),sizeof(double));
  s.read((char*)(&flow),sizeof(double));
  s.read((char*)(&capacity),sizeof(double));
}

