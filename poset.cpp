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

extern Random RND;

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

bool Poset::remove_order(int u,int v)
{
  boost::unordered_map<std::pair<int,int>,bool>::const_iterator qt = order.find(std::pair<int,int>(v,u));
  if (qt == order.end()) return false;

  return true;
}

bool Poset::invert_order(int u,int v)
{
  boost::unordered_map<std::pair<int,int>,bool>::const_iterator qt = order.find(std::pair<int,int>(v,u));
  if (qt == order.end()) return false;

  return true;
}

RELATION Poset::get_order(int u,int v) const
{
  if (u == v) return BEFORE;
  boost::unordered_map<std::pair<int,int>,bool>::const_iterator qt;
  qt = order.find(std::pair<int,int>(u,v));
  if (qt != order.end()) return BEFORE;
  qt = order.find(std::pair<int,int>(v,u));
  if (qt == order.end()) return INCOMPARABLE;
  return AFTER;
}

bool Poset::set_order(int u,int v)
{
  if (u == v) return false;
  // Check to make sure that we don't already have v < u in this poset
  boost::unordered_map<std::pair<int,int>,bool>::const_iterator qt = order.find(std::pair<int,int>(v,u));
  if (qt != order.end()) return false;
  std::pair<int,int> pr(u,v);
  qt = order.find(pr);
  if (qt == order.end()) {
    order[pr] = true;
    // Keep the ordering transitive!
    for(int i=0; i<N; ++i) {
      if (i == u || i == v) continue;
      if (get_order(v,i) == BEFORE) set_order(u,i);
      if (get_order(i,u) == BEFORE) set_order(i,v);
    }
    return true;
  }
  return false;
}

int Poset::chain_number(int length) const
{
  // Compute the number of chains of a given length in this poset 
  int i,j,nchain = 0;
  std::pair<int,int> pr;
  boost::unordered_map<std::pair<int,int>,bool>::const_iterator qt;

  if (length == 2) {
    // The easiest case
    for(i=0; i<N; ++i) {
      pr.first = i;
      for(j=0; j<N; ++j) {
        if (i == j) continue;
        pr.second = j;
        qt = order.find(pr);
        if (qt != order.end()) nchain++;
      }
    }
  }
  else if (length == 3) {
    // Fairly easy as well
    for(i=0; i<N; ++i) {
      pr.first = i;
      for(j=0; j<N; ++j) {
        if (i == j) continue;
        pr.second = j;
        qt = order.find(pr);
        if (qt == order.end()) continue;
        nchain += width(i,j);
      }
    }
  }
  else {
    // The general case, rather complicated
    int k,l;
    for(i=0; i<N; ++i) {
      pr.first = i;
      for(j=0; j<N; ++j) {
        if (i == j) continue;
        pr.second = j;
        qt = order.find(pr);
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
  std::pair<int,int> pr;
  boost::unordered_map<std::pair<int,int>,bool>::const_iterator qt;
  for(i=0; i<N; ++i) {
    if (i == u || i == v) continue;
    pr.first = u; pr.second = i;
    qt = order.find(pr);
    if (qt == order.end()) continue;
    pr.first = i; pr.second = v;
    qt = order.find(pr);
    if (qt == order.end()) continue;
    w++;
  }
  return w;
}

bool Poset::covered(int u,int v) const
{
  // A method to determine if u is covered by v
  if (u == v) return false;
  boost::unordered_map<std::pair<int,int>,bool>::const_iterator qt = order.find(std::pair<int,int>(u,v));
  if (qt == order.end()) return false;    
  return (width(u,v) == 0) ? true : false;
}

bool Poset::consistent() const
{
  // We need to make sure the order relation satisfies the axioms of a 
  // partial order, i.e. reflexive, anti-symmetric and transitive. 
  int i,j,k;
  boost::unordered_map<std::pair<int,int>,bool>::const_iterator qt;

  for(i=0; i<N; ++i) {
    // Find every vertex that is after this one and make sure that all the 
    // vertices after them are also after "i"
    for(j=0; j<N; ++j) {
      if (i == j) continue;
      qt = order.find(std::pair<int,int>(i,j));
      if (qt == order.end()) continue;
      // If i < j, then we cannot have that j < i
      qt = order.find(std::pair<int,int>(j,i));
      if (qt != order.end()) return false;
      for(k=0; k<N; ++k) {
        if (k == j || k == i) continue;
        if (get_order(j,k) == BEFORE) {
          if (get_order(i,k) != BEFORE) return false;
        }
      }
    }
  }
  return true;   
}

void Poset::write_incastrature(const std::string& filename) const
{
  // A method that generates the Hasse diagram corresponding to the
  // complex's simplicial structure, with the diagram stored in the
  // PDF "hasse.pdf"; the method assumes that the Graphviz library
  // has been installed on the system.
  int i,j;

  std::ofstream s(filename.c_str(),std::ios::trunc);

  s << "digraph G {" << std::endl;

  for(i=0; i<N; ++i) {
    for(j=0; j<N; ++j) {
      if (i == j) continue;
      if (get_order(i,j) == AFTER) {
        s << "  \"" << 1+i << "\" -> \"" << 1+j << "\";" << std::endl;
      }
    }
  }
  s << "}" << std::endl;
  s.close();
}

void Poset::construct_order(double lambda)
{
  // A method to impose a random order on the poset 
  int u,v,n = 0;
  double percent;
  boost::unordered_map<std::pair<int,int>,bool>::const_iterator qt;
  const int M = (N*(N-1))/2;

  do {
    u = RND.irandom(N);
    v = RND.irandom(N);
    if (set_order(u,v)) {
      assert(consistent());
      n++;
      percent = double(order.size())/double(M);
    }
  } while(percent < lambda);
}

