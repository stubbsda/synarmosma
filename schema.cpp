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

#include "schema.h"

using namespace SYNARMOSMA;

Schema::Schema()
{
  
}

Schema::Schema(unsigned int n)
{
  unsigned int i;
  for(i=0; i<n; ++i) {
    add_vertex();
  }
}

Schema::Schema(const Schema& source)
{
  vertices = source.vertices;
}

Schema& Schema::operator =(const Schema& source)
{
  if (this == &source) return *this;
  vertices = source.vertices;
  return *this;
}

Schema::~Schema()
{

}

void Schema::clear()
{
  vertices.clear();
}

void Schema::serialize(std::ofstream& s) const
{
  unsigned int i,j,n = get_order();
  bool active;
  std::set<unsigned int> N;
  std::set<unsigned int>::const_iterator it;

  s.write((char*)(&n),sizeof(int));
  for(i=0; i<n; ++i) {
    active = vertices[i].first;    
    s.write((char*)(&active),sizeof(bool));
    N = vertices[i].second;
    j = (signed) N.size();
    s.write((char*)(&j),sizeof(int));
    for(it=N.begin(); it!=N.end(); ++it) {
      j = *it;
      s.write((char*)(&j),sizeof(int));
    }
  }
}

void Schema::deserialize(std::ifstream& s)
{
  unsigned int i,j,k,n,m;
  bool active;
  std::set<unsigned int> N;

  clear();

  s.read((char*)(&n),sizeof(int));
  for(i=0; i<n; ++i) {
    s.read((char*)(&active),sizeof(bool));
    s.read((char*)(&m),sizeof(int));
    for(j=0; j<m; ++j) {
      s.read((char*)(&k),sizeof(int));
      N.insert(k);
    }
    vertices.push_back(std::pair<bool,std::set<unsigned int> >(active,N));
    N.clear();
  }
}

bool Schema::connected(unsigned int n,unsigned int m) const
{
  // A method to check if the vertices n and m share a direct 
  // connection
  if (n == m) return false;
  const int order = get_order();
  assert(n < order);
  assert(m < order);

  if (vertices[n].second.count(m) == 0) {
    // This edge doesn't exist...
    return false;
  }
  return true;
}

unsigned int Schema::add_vertex()
{
  unsigned int nv = vertices.size();
  std::set<int> empty;
  vertices.push_back(std::pair<bool,std::set<int> >(true,empty));
  return nv;
}

bool Schema::drop_vertex(unsigned int n)
{
  const unsigned int order = get_order();
  assert(n < order);
  if (!vertices[n].first) return false;
  unsigned int m;
  std::set<unsigned int> N = vertices[n].second;
  std::set<unsigned int>::const_iterator it,jt;

  vertices[n].first = false;
  vertices[n].second.clear();
  for(it=N.begin(); it!=N.end(); ++it) {
    m = *it;
    jt = vertices[m].second.find(n);
    vertices[m].second.erase(*jt);
  }
  return true;
}

bool Schema::drop_edge(unsigned int n,unsigned int m)
{
  if (n == m) return false;
  const unsigned int order = get_order();
  assert(n < order);
  assert(m < order);

  std::set<unsigned int>::const_iterator it;

  it = vertices[n].second.find(m);
  if (it == vertices[n].second.end()) {
    // This edge doesn't exist...
    return false;
  }
  vertices[n].second.erase(*it);
  it = vertices[m].second.find(n);
  vertices[m].second.erase(*it);
  return true;
}

bool Schema::add_edge(unsigned int n,unsigned int m)
{
  if (n == m) return false;
  const unsigned int order = get_order();
  assert(n < order);
  assert(m < order);

  if (vertices[n].second.count(m) == 0) {
    // This edge doesn't already exist, so add it...
    vertices[n].second.insert(m);
    vertices[m].second.insert(n);
    return true;
  }
  return false;
}

bool Schema::positive_valence() const
{
  // This method just checks if there are any active vertices with no connections,
  // in which case it returns false, true otherwise
  unsigned int i,order = get_order();
  for(i=0; i<order; ++i) {
    if (vertices[i].first && vertices[i].second.empty()) return false;
  }
  return true;
}

bool Schema::consistent() const
{
  // Another utility routine, useful in debugging, that checks for various
  // pathological conditions in the neighbour table, such as vertices that have
  // a valence less than or equal to zero, or claim to be connected to themselves,
  // or have two bonds to the same vertex, or have (locally) inconsistent bonds,
  // or claim to be bonded to nonexistent vertices. If any of these conditions
  // arise, then a relatively verbose error message is printed out, and the method
  // returns false.
  int i,in1;
  std::set<int>::const_iterator it;

  // Loop through all vertices
  for(i=0; i<nvertex; ++i) {
    // Check if the valence is weird
    if (neighbours[i].size() == 0) return false;
    // Now loop through all the neighbours
    for(it=neighbours[i].begin(); it!=neighbours[i].end(); ++it) {
      in1 = *it;
      // Check for self-bonding
      if (in1 == i) return false;
      // Check for neighbour values that are negative or too large
      if (in1 >= nvertex) return false;
      // If this neighbour vertex is local, then check to see if it too lists
      // i as a neighbour
      if (neighbours[in1].count(i) != 1) {
        std::cerr << "Schema inconsistent at " << i << " and " << in1 << std::endl;
        std::exit(1);
      }
    }
  }
  return true;
}

void Schema::components(std::vector<int>& csize,std::vector<int>& celements) const
{
  int i,start,nc,noe;
  std::vector<int> component;
  std::set<int> change,nchange;
  std::set<int>::const_iterator it1,it2;

  nc = 0;
  for(i=0; i<nvertex; ++i) {
    component.push_back(-1);
  }

  do {
    start = -1;
    for(i=0; i<nvertex; ++i) {
      if (component[i] == -1) {
        start = i;
        break;
      }
    }
    if (start == -1) break;
    change.insert(start);
    noe = 0;
    do {
      for(it1=change.begin(); it1!=change.end(); it1++) {
        component[*it1] = nc;
        celements.push_back(*it1);
      }
      noe += change.size();
      for(it1=change.begin(); it1!=change.end(); it1++) {
        i = *it1;
        for(it2=neighbours[i].begin(); it2!=neighbours[i].end(); it2++) {
          if (component[*it2] < 0) nchange.insert(*it2);
        }
      }
      if (nchange.empty()) break;
      change = nchange;
      nchange.clear();
    } while(true);
    csize.push_back(noe);
    nc++;
    change.clear();
  } while(true);
}

bool Schema::connected() const
{
  int i,n;
  std::vector<int> ubiquity;
  std::set<int> change,nchange;
  std::set<int>::const_iterator it,jt;

  for(i=0; i<nvertex; ++i) {
    ubiquity.push_back(0);
  }
  ubiquity[0] = 1;
  change.insert(0);

  do {
    for(it=change.begin(); it!=change.end(); ++it) {
      i = *it;
      for(jt=neighbours[i].begin(); jt!=neighbours[i].end(); ++jt) {
        n = *jt;
        if (ubiquity[n] == 1) continue;
        nchange.insert(n);
      }
    }
    if (nchange.empty()) break;
    for(it=nchange.begin(); it!=nchange.end(); ++it) {
      n = *it;
      ubiquity[n] = 1;
    }
    change = nchange;
    nchange.clear();
  } while(true);
  for(i=0; i<nvertex; ++i) {
    if (ubiquity[i] == 0) return false;
  }
  return true;
}

int Schema::component_analysis(std::vector<int>& component) const
{
  int i,n,m,start,nc = 0;
  std::set<int> change,nchange;
  std::set<int>::const_iterator it,jt;

  component.clear();
  for(i=0; i<nvertex; ++i) {
    component.push_back(-1);
  }

  do {
    start = -1;
    for(i=0; i<nvertex; ++i) {
      if (component[i] == -1) {
        start = i;
        break;
      }
    }
    if (start == -1) break;
    change.insert(start);
    do {
      for(it=change.begin(); it!=change.end(); ++it) {
        n = *it;
        component[n] = nc;
      }
      for(it=change.begin(); it!=change.end(); ++it) {
        n = *it;
        for(jt=neighbours[n].begin(); jt!=neighbours[n].end(); ++jt) {
          m = *jt;
          if (component[m] < 0) nchange.insert(m);
        }
      }
      if (nchange.empty()) break;
      change = nchange;
      nchange.clear();
    } while(true);
    nc++;
    change.clear();
  } while(true);
  return nc;
}

int Schema::spanning_tree(std::vector<int>& tree_edges) const
{
  int ntree,n,m;
  std::set<int> current,vertices,next;
  std::set<int>::const_iterator it,jt,kt,lt;

  // A sanity check...
  assert(connected());

  current.insert(0);

  do {
    for(it=current.begin(); it!=current.end(); ++it) {
      vertices.insert(*it);
    }
    for(it=current.begin(); it!=current.end(); ++it) {
      n = *it;
      for(jt=neighbours[n].begin(); jt!=neighbours[n].end(); ++jt) {
        m = *jt;
        kt = std::find(vertices.begin(),vertices.end(),m);
        if (kt != vertices.end()) continue;
        lt = std::find(next.begin(),next.end(),m);
        if (lt != next.end()) continue;
        next.insert(m);
        // Now add this edge to the spanning tree...
        if (n < m) {
          tree_edges.push_back(n); tree_edges.push_back(m);
        }
        else {
          tree_edges.push_back(m); tree_edges.push_back(n);
        }
      }
    }
    if (next.empty()) break;
    current = next;
    next.clear();
  } while(true);
  ntree = (signed) tree_edges.size();
  assert(nvertex == (1+ntree/2));
  return ntree;
}



