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
}

int Vertex::serialize(std::ofstream& s) const
{
  int i,n,count = 0;
  std::set<int>::const_iterator it;

  s.write((char*)(&incept),sizeof(int)); count += sizeof(int);
  s.write((char*)(&topological_dimension),sizeof(int)); count += sizeof(int);
#ifdef DISCRETE
  s.write((char*)(&energy),sizeof(UINT64)); count += sizeof(UINT64);
#else
  s.write((char*)(&energy),sizeof(double)); count += sizeof(double);
#endif
  n = (signed) anterior.size();
  s.write((char*)(&n),sizeof(int)); count += sizeof(int);
  for(it=anterior.begin(); it!=anterior.end(); ++it) {
    i = *it;
    s.write((char*)(&i),sizeof(int)); count += sizeof(int);
  }
  n = (signed) posterior.size();
  s.write((char*)(&n),sizeof(int)); count += sizeof(int);
  for(it=posterior.begin(); it!=posterior.end(); ++it) {
    i = *it;
    s.write((char*)(&i),sizeof(int)); count += sizeof(int);
  }
  n = (signed) neighbours.size();
  s.write((char*)(&n),sizeof(int)); count += sizeof(int);
  for(it=neighbours.begin(); it!=neighbours.end(); ++it) {
    i = *it;
    s.write((char*)(&i),sizeof(int)); count += sizeof(int);
  }
  n = (signed) entourage.size();
  s.write((char*)(&n),sizeof(int)); count += sizeof(int);
  for(it=entourage.begin(); it!=entourage.end(); ++it) {
    i = *it;
    s.write((char*)(&i),sizeof(int)); count += sizeof(int);
  }
  return count;
}

int Vertex::deserialize(std::ifstream& s)
{
  int i,j,n,count = 0;

  clear();

  s.read((char*)(&incept),sizeof(int)); count += sizeof(int);
  s.read((char*)(&topological_dimension),sizeof(int)); count += sizeof(int);
#ifdef DISCRETE
  s.read((char*)(&energy),sizeof(UINT64)); count += sizeof(UINT64);
#else
  s.read((char*)(&energy),sizeof(double)); count += sizeof(double);
#endif
  s.read((char*)(&n),sizeof(int)); count += sizeof(int);
  for(i=0; i<n; ++i) {
    s.read((char*)(&j),sizeof(int)); count += sizeof(int);
    anterior.insert(j);
  }
  s.read((char*)(&n),sizeof(int)); count += sizeof(int);
  for(i=0; i<n; ++i) {
    s.read((char*)(&j),sizeof(int)); count += sizeof(int);
    posterior.insert(j);
  }
  s.read((char*)(&n),sizeof(int)); count += sizeof(int);
  for(i=0; i<n; ++i) {
    s.read((char*)(&j),sizeof(int)); count += sizeof(int);
    neighbours.insert(j);
  }
  s.read((char*)(&n),sizeof(int)); count += sizeof(int);
  for(i=0; i<n; ++i) {
    s.read((char*)(&j),sizeof(int)); count += sizeof(int);
    entourage.insert(j);
  }
  return count;
}









