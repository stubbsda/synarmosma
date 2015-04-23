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

#include "poset.h"

using namespace SYNARMOSMA;

Poset::Poset()
{
  N = 0;
}

Poset::Poset(int n)
{
  // We begin by assuming that every pair of points is 
  // incomparable, so the hash map is empty
  N = n;
}

Poset::~Poset()
{

}

void Poset::clear()
{
  N = 0;
  order.clear();
}

void Poset::add_vertex()
{
  N += 1;
}

void Poset::set_order(int u,int v)
{
  order[u,v] = true;
}

int Poset::chain_number(int length) const
{
  // Compute the number of chains of a given length in this poset 
  int i,j,nchain = 0;
  boost::unordered_map<int,int,bool>::const_iterator qt;

  if (length == 2) {
    // The easiest case
    for(i=0; i<N; ++i) {
      for(j=1+i; j<N; ++j) {
        qt = order.find(i,j);
        if (qt != order.end()) nchain++;
      }
    }
  }
  else if (length == 3) {
    // Fairly easy as well
    for(i=0; i<N; ++i) {
      for(j=1+i; j<N; ++j) {
        qt = order.find(i,j);
        if (qt == order.end()) continue;
        nchain += width(i,j);
      }
    }
  }
  else {
    // The general case, rather complicated
    int k,l;
    for(i=0; i<N; ++i) {
      for(j=1+i; j<N; ++j) {
        qt = order.find(i,j);
        if (qt == order.end()) continue;
        // See how many chains of length l between the elements i and j
        // can be built
        l = width(i,j);
        for(k=0; k<length-2; ++k) {

        }
      }
    }
  }
  return nchain;
}

int Poset::width(int u,int v) const
{
  assert(u != v);
  int i,w = 0;
  boost::unordered_map<int,int,bool>::const_iterator qt;
  for(i=0; i<N; ++i) {
    if (i == u || i == v) continue;
    qt = order.find(u,i);
    if (qt == order.end()) continue;
    qt = order.find(i,v);
    if (qt == order.end()) continue;
    w++;
  }
  return w;
}

bool Poset::covered(int u,int v) const
{
  // A method to determine if u is covered by v
  if (u == v) return false;
  boost::unordered_map<int,int,bool>::const_iterator qt = order.find(u,v);
  if (qt == order.end()) return false;    
  return (width(u,v) == 0) ? true : false;
}

RELATION Poset::get_relation(int u,int v) const
{
  if (u == v) return BEFORE;
  boost::unordered_map<int,int,bool>::const_iterator qt;

  qt = order.find(u,v);
  if (qt != order.end()) return BEFORE;
  qt = order.find(v,u);
  if (qt == order.end()) return INCOMPARABLE;
  return AFTER;
}

bool Poset::consistent() const
{
  // We need to make sure the order relation satisfies the axioms of a 
  // partial order, i.e. reflexive, anti-symmetric and transitive. The 
  // first two we obtain automatically, let's check the last one
  int i,j,k;
  boost::unordered_map<int,int,bool>::const_iterator qt;

  for(i=0; i<N; ++i) {
    // Find every vertex that is after this one and make sure that all the 
    // vertices after them are also after "i"
    for(j=1+i; j<N; ++j) {
      qt = order.find(i,j);
      if (qt == order.end()) continue;
      for(k=0; k<N; ++k) {
        if (k == j || k == i) continue;
        if (get_relation(j,k) == BEFORE) {
          if (get_relation(i,k) != BEFORE) return false;
        }
      }
    }
  }
  return true;   
}
