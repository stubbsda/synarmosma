#include "lattice.h"

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

