#include "edge.h"

Edge::Edge()
{

}

Edge::~Edge()
{

}

Edge::Edge(const Edge& source)
{
  length = source.length;
  direction = source.direction;
  flow = source.flow;
  capacity = source.capacity;
  colour = source.colour;
  nodes[0] = source.nodes[0];
  nodes[1] = source.nodes[1];
  name = source.name;
  cyclic = source.cyclic;
}

Edge& Edge::operator =(const Edge& source)
{
  if (this == &source) return *this;

  length = source.length;
  direction = source.direction;
  flow = source.flow;
  capacity = source.capacity;
  colour = source.colour;
  nodes[0] = source.nodes[0];
  nodes[1] = source.nodes[1];
  name = source.name;
  cyclic = source.cyclic;

  return *this;
}

void Edge::clear()
{
  length = 0.0;
  nodes[0] = 0;
  nodes[1] = 0;
  capacity = 1.0;
  flow = 0.0;
  colour = 2;
  name = "";
  direction = 0;
  cyclic = false;
}

