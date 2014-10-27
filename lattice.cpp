#include "lattice.h"

using namespace SYNARMOSMA;

Lattice::Lattice()
{
  N = 0;
}

Lattice::Lattice(int n)
{
  int i,j;

  N = n;
  // We begin by assuming that every pair of points is 
  // incomparable
  for(i=0; i<N; ++i) {
    for(j=1+i; j<N; ++j) {
      order[make_key(i,j)] = INCOMPARABLE;
    }
  }
}

Lattice::~Lattice()
{

}

void Lattice::clear()
{
  N = 0;
  order.clear();
}

void Lattice::add_vertex()
{
  N += 1;
  for(int i=0; i<N-1; ++i) {
    order[make_key(i,N-1)] = INCOMPARABLE;
  }
}

RELATION Lattice::get_relation(int u,int v) const
{
  if (u == v) return BEFORE; 
  boost::unordered_map<std::string,RELATION>::const_iterator qt = order.find(make_key(u,v));
  if (qt->second == INCOMPARABLE) return INCOMPARABLE;
  if (u < v) {
    return qt->second;
  }
  else {
    RELATION output = (qt->second == BEFORE) ? AFTER : BEFORE;
    return output;
  }
}

bool Lattice::consistent() const
{
  // We need to make sure the order relation satisfies the axioms of a 
  // partial order, i.e. reflexive, anti-symmetric and transitive. The 
  // first two we obtain automatically, let's check the last one
  int i,j,k;
  boost::unordered_map<std::string,RELATION>::const_iterator qt;

  for(i=0; i<N; ++i) {
    // Find every vertex that is after this one and make sure that all the 
    // vertices after them are also after "i"
    for(j=1+i; j<N; ++j) {
      qt = order.find(make_key(i,j));
      if (qt->second == BEFORE) {
        for(k=0; k<N; ++k) {
          if (k == j || k == i) continue;
          if (get_relation(j,k) == BEFORE) {
            if (get_relation(i,k) != BEFORE) return false;
          }
        }
      }
    }
  }
  return true;   
}
