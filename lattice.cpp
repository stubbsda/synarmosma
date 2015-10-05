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

extern Random RND;

Lattice::Lattice() : Poset()
{
  atomic = false;
  null = 0;
  unity = 0;
  initialize();
}

Lattice::Lattice(unsigned int n) : Poset(n)
{
  atomic = false;
  null = 0;
  unity = 0;
  initialize();
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
  N = 0;
  order.clear();
  atomic = false;
  null = 0;
  unity = 0;
  atoms.clear();
}

void Lattice::serialize(std::ofstream& s) const
{
  unsigned int i,j;
  std::set<unsigned int>::const_iterator it;
  RELATION rho;

  s.write((char*)(&N),sizeof(int));
  s.write((char*)(&atomic),sizeof(bool));
  s.write((char*)(&null),sizeof(int));
  s.write((char*)(&unity),sizeof(int));
  for(i=0; i<N; ++i) {
    for(j=1+i; j<N; ++j) {
      rho = get_order(i,j);
      s.write((char*)(&rho),sizeof(int));
    }
  }
  j = atoms.size();
  s.write((char*)(&j),sizeof(int));
  for(it=atoms.begin(); it!=atoms.end(); ++it) {
    i = *it;
    s.write((char*)(&i),sizeof(int));
  }
}

void Lattice::deserialize(std::ifstream& s)
{
  unsigned int i,j,k;
  RELATION rho;

  clear();

  s.read((char*)(&N),sizeof(int));
  s.read((char*)(&atomic),sizeof(bool));
  s.read((char*)(&null),sizeof(int));
  s.read((char*)(&unity),sizeof(int));
  for(i=0; i<N; ++i) {
    for(j=1+i; j<N; ++j) {
      s.read((char*)(&rho),sizeof(int));
      if (rho == BEFORE) {
        set_order(i,j);
      }
      else if (rho == AFTER) {
        set_order(j,i);
      }
    }
  }
  s.read((char*)(&j),sizeof(int));
  for(i=0; i<j; ++i) {
    s.read((char*)(&k),sizeof(int));
    atoms.insert(k);
  }
  assert(consistent());
}


void Lattice::initialize()
{
  if (N < 2) return; 
  unsigned int i,j,n1,n2,delta = 2*N*(N-1),ndelta;

  do {
    n1 = RND.irandom(N);
    n2 = RND.irandom(N);
    if (n1 == n2) continue;
    if (get_order(n1,n2) != DISPARATE) continue;
    set_order(n1,n2);
    
    ndelta = 0;
    for(i=0; i<N; ++i) {
      for(j=0; j<N; ++j) {
        if (i == j) continue;
        if (meet(i,j) == N) ndelta++;
        if (join(i,j) == N) ndelta++;
      }
    }
    if (ndelta > delta) {
      unset_order(n1,n2);
      continue;
    }
    delta = ndelta;
  } while(delta > 0);
  assert(consistent());
  compute_bounds();
  compute_atoms();
}

bool Lattice::consistent() const
{
  // In a lattice every pair of elements must have a meet and join, so...
  unsigned int i,j;

  for(i=0; i<N; ++i) {
    for(j=0; j<N; ++j) {
      if (i == j) continue;
      if (meet(i,j) == N) return false;
      if (join(i,j) == N) return false;
    }
  }
  return true;
}

void Lattice::compute_bounds()
{
  // This method computes the null and identity elements for this lattice using a 
  // brute force technique
  unsigned int i,j;
  bool found,nfound;
  
  // We begin by finding the 0 element, i.e. the element of the lattice such that 0 <= x for all x
  nfound = false;
  for(i=0; i<N; ++i) {
    found = true;
    for(j=0; j<N; ++j) {
      if (i == j) continue;
      if (get_order(i,j) != BEFORE) {
        found = false;
        break;
      }
    }
    if (found) {
      null = i;
      nfound = true;
      break;
    }
  }
  if (!nfound) {    
    null = join(0,1);
    for(i=0; i<N-1; ++i) {
      null = join(null,i+1);
    }
    N += 1;
  }

  // Now the 1 element, i.e. the element of the lattice such that x <= 1 for all x
  nfound = false;
  for(i=0; i<N; ++i) {
    found = true;
    for(j=0; j<N; ++j) {
      if (i == j) continue;
      if (get_order(i,j) != AFTER) {
        found = false;
        break;
      }
    }
    if (found) {
      unity = i;
      nfound = true;
      break;
    }
  }
  if (!nfound) {
    unity = meet(0,1);
    for(i=0; i<N-1; ++i) {
      unity = meet(unity,i+1);
    }
    N += 1;
  }   
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
  unsigned int i,j,rvalue = N;
  bool max;
  std::set<unsigned int> candidates;
  std::set<unsigned int>::const_iterator it,jt;
  RELATION rho;
  
  for(i=0; i<N; ++i) {
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
      if (rho == BEFORE) {
        max = false;
        break;
      }
    }
    if (max) return i;
  }
  return rvalue;
}

unsigned int Lattice::join(unsigned int x,unsigned int y) const
{
  // Find the element w in L such that x <= w and y <= w, while any other 
  // z in L satisfying these relations is such that w <= z.
  unsigned int i,j,rvalue = N;
  bool max;
  std::set<unsigned int> candidates;
  std::set<unsigned int>::const_iterator it,jt;
  RELATION rho;
  
  for(i=0; i<N; ++i) {
    rho = get_order(i,x);
    if (rho != AFTER && i != x) continue;
    rho = get_order(i,y);
    if (rho != AFTER && i != y) continue;
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
      if (rho == AFTER) {
        max = false;
        break;
      }
    }
    if (max) return i;
  }
  return rvalue;
}

