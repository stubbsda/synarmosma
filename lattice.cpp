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
  atomic = false;
  null = 0;
  unity = 0;
}

Lattice::Lattice(unsigned int n) : Poset(n)
{
  atomic = false;
  null = 0;
  unity = 0;
  compute_bounds();
  compute_atoms();
}

Lattice::~Lattice()
{

}

Lattice::Lattice(const Lattice& source)
{
  N = source.N;
  order = source.order;
  atomic = source.atomic;
  atoms = source.atoms;
  null = source.null;
  unity = source.unity;
}

Lattice& Lattice::operator =(const Lattice& source) 
{
  if (this == &source) return *this;
  N = source.N;
  order = source.order;
  atomic = source.atomic;
  atoms = source.atoms;
  null = source.null;
  unity = source.unity;
  return *this;
}

void Lattice::clear()
{
  Poset::clear();
  atomic = false;
  null = 0;
  unity = 0;
  atoms.clear();
}

bool Lattice::consistent() const
{
  // In a lattice every pair of elements must have a meet and join, so...
  unsigned int i,j,n;

  for(i=0; i<N; ++i) {
    for(j=0; j<N; ++j) {
      if (i == j) continue;
      n = meet(i,j);
      n = join(i,j);
    }
  }
  return true;
}

void Lattice::compute_bounds()
{
  // This method computes the null and identity elements for this lattice
}

void Lattice::compute_atoms()
{
  unsigned int i,j;
  bool found;
  std::set<unsigned int>::const_iterator it;

  atoms.clear();
  for(i=0; i<N; ++i) {
    if (i == null) continue;
    if (covered(null,i)) atoms.insert(i);
  }
  
  atomic = true;
  for(i=0; i<N; ++i) {
    if (i == null) continue;
    if (atoms.count(i) > 0) continue;
    found = false;
    for(it=atoms.begin(); it!=atoms.end(); ++it) {
      j = *it;
      if (get_order(j,i) == BEFORE) {
        found = true;
        break;
      }
    }
    if (!found) {
      atomic = false;
      break;
    }
  }
}

unsigned int Lattice::meet(unsigned int x,unsigned int y) const
{
  // Find the element w in L such that w <= x and w <= y, while any other 
  // z in L satisfying these relations is such that z <= w.
  unsigned int i,j;
  bool max;
  std::set<unsigned int> candidates;
  std::set<unsigned int>::const_iterator it,jt;
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
  std::cerr << "Lattice pair which doesn't have a meet, exiting!" << std::endl;
  std::exit(1);
}

unsigned int Lattice::join(unsigned int x,unsigned int y) const
{
  // Find the element w in L such that x <= w and y <= w, while any other 
  // z in L satisfying these relations is such that w <= z.
  unsigned int i,j;
  bool max;
  std::set<unsigned int> candidates;
  std::set<unsigned int>::const_iterator it,jt;
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
  std::cerr << "Lattice pair which doesn't have a meet, exiting!" << std::endl;
  std::exit(1);
}

