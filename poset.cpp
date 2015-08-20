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

Poset::Poset(const Poset& source)
{
  N = source.N;
  order = source.order;
}

Poset& Poset::operator =(const Poset& source) 
{
  if (this == &source) return *this;
  N = source.N;
  order = source.order;
  return *this;
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

bool Poset::invert_order(int u,int v)
{
  if (u == v) return false;
  RELATION rho = get_order(u,v);
  if (rho == INCOMPARABLE) return false;
  if (rho == BEFORE) {
    // Inverting order of u < v
    set_order(v,u);
  }
  else {
    // Inverting order of v < u
    set_order(u,v);
  }
  return true;
}

bool Poset::unset_order(int u,int v)
{
  if (u == v) return false;
  RELATION rho = get_order(u,v);
  if (rho == INCOMPARABLE) return false;
  if (rho == BEFORE) {
    // Removing order of u < v
    order.erase(std::pair<int,int>(u,v));
    for(int i=0; i<N; ++i) {
      if (i == u || i == v) continue;
      if (get_order(u,i) == BEFORE && get_order(i,v) == BEFORE) unset_order(i,v);
    }
  }
  else {
    // Removing order of v < u
    order.erase(std::pair<int,int>(v,u));
    for(int i=0; i<N; ++i) {
      if (i == u || i == v) continue;
      if (get_order(v,i) == BEFORE && get_order(i,u) == BEFORE) unset_order(i,u);
    }
  }
  return true;
}

bool Poset::set_order(int u,int v)
{
  if (u == v) return false;
  // Check to see if we already have u < v or v < u in this poset
  RELATION rho = get_order(u,v);
  if (rho == BEFORE) return false;
  if (rho == AFTER) order.erase(std::pair<int,int>(v,u));
  order[std::pair<int,int>(u,v)] = true;
  // Keep the ordering transitive!
  for(int i=0; i<N; ++i) {
    if (i == u || i == v) continue;
    if (get_order(v,i) == BEFORE) set_order(u,i);
    if (get_order(i,u) == BEFORE) set_order(i,v);
  }
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

int Poset::build_chain(std::vector<int>& chain,int length) const
{
  int l = (signed) chain.size(),output = 0;
  if (l == length) {
    output = 1;
  }
  else {
    std::vector<int> nchain = chain;

    for(int i=0; i<N; ++i) {
      if (i == chain[l]) continue;
      if (get_order(chain[l],i) == BEFORE) {
        nchain.push_back(i);
        output += build_chain(nchain,length);
        nchain = chain;
      }
    }
  }
  return output;
}

int Poset::chain_number(int length) const
{
  // Compute the number of chains of a given length in this poset 
  assert(length >= 0);
  if (length == 0 || length == 1) return 0;
  int i,j,nchain = 0;
  std::set<int> sigma;

  if (length == 2) {
    // The easiest case
    for(i=0; i<N; ++i) {
      for(j=0; j<N; ++j) {
        if (i == j) continue;
        if (get_order(i,j) == BEFORE) nchain++;
      }
    }
  }
  else if (length == 3) {
    // Fairly easy as well
    for(i=0; i<N; ++i) {
      for(j=0; j<N; ++j) {
        if (i == j) continue;
        if (get_order(i,j) != BEFORE) continue;
        compute_width(i,j,sigma);
        nchain += (signed) sigma.size();
      }
    }
  }
  else {
    // The general case, rather complicated so we use recursion to do it
    std::vector<int> chain;
   
    for(i=0; i<N; ++i) {
      chain.push_back(i);
      nchain += build_chain(chain,length);
      chain.clear();
    }   
  }
  return nchain;
}

void Poset::compute_width(int u,int v,std::set<int>& slice) const
{
  assert(u != v);

  slice.clear();

  for(int i=0; i<N; ++i) {
    if (i == u || i == v) continue;
    if (get_order(u,i) != BEFORE) continue;
    if (get_order(i,v) != BEFORE) continue;
    slice.insert(i);
  }
}

bool Poset::covered(int u,int v) const
{
  // A method to determine if u is covered by v
  if (u == v) return false;
  if (get_order(u,v) != BEFORE) return false;
  std::set<int> sigma;
  compute_width(u,v,sigma);
  return sigma.empty(); 
}

bool Poset::sink(int n) const
{
  // This is an element whose posteriority is null
  std::set<int> S;
  compute_posteriority(n,S);
  return S.empty();
}

bool Poset::source(int n) const
{
  // This is an element whose anteriority is null
  std::set<int> S;
  compute_anteriority(n,S);
  return S.empty();
}

void Poset::compute_anteriority(int n,std::set<int>& output) const
{
  output.clear();
  for(int i=0; i<N; ++i) {
    if (i == n) continue;
    if (get_order(i,n) == BEFORE) output.insert(i);
  }
}

void Poset::compute_posteriority(int n,std::set<int>& output) const
{
  output.clear();
  for(int i=0; i<N; ++i) {
    if (i == n) continue;
    if (get_order(n,i) == BEFORE) output.insert(i);
  }
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
      if (get_order(i,j) != BEFORE) continue;
      if (get_order(j,i) != AFTER) return false;
      for(k=0; k<N; ++k) {
        if (k == j || k == i) continue;
        // So if i < j and j < k, then it must be that i < k
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
  // First all the elements in the poset...
  for(i=0; i<N; ++i) {
    s << "  \"" << 1+i << "\";" << std::endl;
  }
  // Now the directed edges induced by the poset's ordering...
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

void Poset::construct_ordering(double lambda)
{
  // A method to impose a random order on the poset 
  int u,v,n = 0;
  double percent;
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

