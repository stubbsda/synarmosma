#include "edge.h"

using namespace SYNARMOSMA;

Edge::Edge()
{
  clear();
}

Edge::Edge(int u,int v,double ell,int d)
{
#ifdef DEBUG
  assert(u != v);
#endif
  if (u < v) {
    low = u;
    high = v;
    direction = d;
  }
  else {
    low = v;
    high = u;
    direction = -d;
  }
  capacity = ell;
}

Edge::~Edge()
{

}

Edge::Edge(const Edge& source)
{
  length = source.length;
  flow = source.flow;
  capacity = source.capacity;
  low = source.low;
  high = source.high;
  cyclic = source.cyclic;
  direction = source.direction;
}

Edge& Edge::operator =(const Edge& source)
{
  if (this == &source) return *this;

  length = source.length;
  flow = source.flow;
  capacity = source.capacity;
  low = source.low;
  high = source.high;
  cyclic = source.cyclic;
  direction = source.direction;

  return *this;
}

void Edge::clear()
{
  low = -1;
  high = -1;
  length = 0.0;
  capacity = 0.0;
  flow = 0.0;
  cyclic = false;
  direction = UNDIRECTED;
}

void Edge::serialize(std::ofstream& s) const
{
  s.write((char*)(&low),sizeof(int));
  s.write((char*)(&high),sizeof(int));
  s.write((char*)(&cyclic),sizeof(bool));
  s.write((char*)(&length),sizeof(double));
  s.write((char*)(&flow),sizeof(double));
  s.write((char*)(&capacity),sizeof(double));
  s.write((char*)(&direction),sizeof(int));
}

void Edge::deserialize(std::ifstream& s)
{
  clear();

  s.read((char*)(&low),sizeof(int));
  s.read((char*)(&high),sizeof(int));
  s.read((char*)(&cyclic),sizeof(bool));
  s.read((char*)(&length),sizeof(double));
  s.read((char*)(&flow),sizeof(double));
  s.read((char*)(&capacity),sizeof(double));
  s.read((char*)(&direction),sizeof(int));
}

