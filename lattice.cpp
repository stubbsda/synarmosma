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

#include "lattice.h"

using namespace SYNARMOSMA;

Lattice::Lattice() : Poset()
{

}

Lattice::Lattice(int n) : Poset(n)
{

}

Lattice::~Lattice()
{

}

Lattice::Lattice(const Lattice& source)
{
  N = source.N;
  order = source.order;
}

Lattice& Lattice::operator =(const Lattice& source) 
{
  if (this == &source) return *this;
  N = source.N;
  order = source.order;
  return *this;
}

void Lattice::clear()
{
  Poset::clear();
}

bool Lattice::consistent() const
{
  // In a lattice every pair of elements must have a meet and join, so...
  int i,j;

  for(i=0; i<N; ++i) {
    for(j=0; j<N; ++j) {
      if (i == j) continue;
      if (meet(i,j) == -1) return false;
      if (join(i,j) == -1) return false;
    }
  }
  return true;
}

void Lattice::add_vertex()
{
  Poset::add_vertex();
}

int Lattice::meet(int x,int y) const
{
  // Find the element w in L such that w <= x and w <= y, while any other 
  // z in L satisfying these relations is such that z <= w.
  int i,j;
  bool max;
  std::set<int> candidates;
  std::set<int>::const_iterator it,jt;
  RELATION rho;
  
  for(i=0; i<N; ++i) {
    if (i == x || i == y) continue;
    rho = get_order(i,x);
    if (rho != BEFORE) continue;
    rho = get_order(i,y);
    if (rho != BEFORE) continue;
    candidates.insert(i); 
  }
  // Find which among the candidates is the largest
  for(it=candidates.begin(); it!=candidates.end(); ++it) {
    i = *it;
    max = true;
    for(jt=candidates.begin(); jt!=candidates.end(); ++jt) {
      j = *jt;
      if (i == j) continue;
      rho = get_order(i,j);
      if (rho != AFTER) {
        max = false;
        break;
      }
    }
    if (max) return i;
  }
  return -1;
}

int Lattice::join(int x,int y) const
{
  // Find the element w in L such that x <= w and y <= w, while any other 
  // z in L satisfying these relations is such that w <= z.
  int i,j;
  bool max;
  std::set<int> candidates;
  std::set<int>::const_iterator it,jt;
  RELATION rho;
  
  for(i=0; i<N; ++i) {
    if (i == x || i == y) continue;
    rho = get_order(i,x);
    if (rho != AFTER) continue;
    rho = get_order(i,y);
    if (rho != AFTER) continue;
    candidates.insert(i); 
  }
  // Find which among the candidates is the largest
  for(it=candidates.begin(); it!=candidates.end(); ++it) {
    i = *it;
    max = true;
    for(jt=candidates.begin(); jt!=candidates.end(); ++jt) {
      j = *jt;
      if (i == j) continue;
      rho = get_order(i,j);
      if (rho != BEFORE) {
        max = false;
        break;
      }
    }
    if (max) return i;
  }
  return -1;
}

