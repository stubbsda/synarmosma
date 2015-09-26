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

#include "directed_graph.h"

using namespace SYNARMOSMA;

extern Random RND;

Directed_Graph::Directed_Graph() : Graph() 
{

}

Directed_Graph::Directed_Graph(int n) : Graph(n)
{
  for(int i=0; i<n; ++i) {
    orientation.add_element();
  }
}

Directed_Graph::Directed_Graph(int n,double rho) : Graph()
{
  // Something probabilistic
  int i,u,v;
  double p = 0.0;
  const double MX = double(n*(n-1))/2.0;

  for(i=0; i<n; ++i) {
    add_vertex();
  }

  do {
    u = RND.irandom(nvertex);
    v = RND.irandom(nvertex);
    if (u == v) continue;
    if (add_edge(u,v)) p = double(edges.size())/MX;  
  } while(p < rho);
  assert(connected());
}

Directed_Graph::~Directed_Graph()
{

}

Directed_Graph::Directed_Graph(const Directed_Graph& source)
{
  nvertex = source.nvertex;
  neighbours = source.neighbours;
  edges = source.edges;
  index_table = source.index_table;
  orientation = source.orientation;
}

Directed_Graph& Directed_Graph::operator =(const Directed_Graph& source) 
{
  if (this == &source) return *this;
  nvertex = source.nvertex;
  neighbours = source.neighbours;
  edges = source.edges;
  index_table = source.index_table;
  orientation = source.orientation;
  return *this;
}

void Directed_Graph::clear()
{
  nvertex = 0;
  neighbours.clear();
  edges.clear();
  index_table.clear();
  orientation.clear();
}

bool Directed_Graph::consistent() const
{
  if (!Graph::consistent()) return false;
  return orientation.consistent();
}

int Directed_Graph::directedness() const
{
  // A method to calculate how many of the edges are directed...
  int vx[2],output,null = 0;
  std::vector<Edge>::const_iterator it;

  for(it=edges.begin(); it!=edges.end(); ++it) {
    if (it->active) {
      it->get_nodes(vx);
      if (orientation.get_order(vx[0],vx[1]) == INCOMPARABLE) null++;
    }
  }
  output = size() - null;
  return output;
}

void Directed_Graph::serialize(std::ofstream& s) const
{
  int i,j,n;
  std::set<int>::const_iterator it;

  s.write((char*)(&nvertex),sizeof(int));
  for(i=0; i<nvertex; ++i) {
    j = (signed) neighbours[i].size();
    s.write((char*)(&j),sizeof(int));
    for(it=neighbours[i].begin(); it!=neighbours[i].end(); ++it) {
      j = *it;
      s.write((char*)(&j),sizeof(int));
    }
  }
  n = (signed) edges.size();
  s.write((char*)(&n),sizeof(int));
  for(i=0; i<n; ++i) {
    edges[i].serialize(s);
  }
  orientation.serialize(s);
}

void Directed_Graph::deserialize(std::ifstream& s)
{
  int i,j,k,n;
  std::set<int> S;
  Edge q;

  clear();

  s.read((char*)(&nvertex),sizeof(int));
  for(i=0; i<nvertex; ++i) {
    s.read((char*)(&n),sizeof(int));
    for(j=0; j<n; ++j) {
      s.read((char*)(&k),sizeof(int));
      S.insert(k);
    }
    neighbours.push_back(S);
    S.clear();
  }
  s.read((char*)(&n),sizeof(int));
  for(i=0; i<n; ++i) {
    q.deserialize(s);
    edges.push_back(q);
  }
  orientation.deserialize(s);
}

int Directed_Graph::add_vertex()
{
  int n = Graph::add_vertex();
  orientation.add_element();
  return n;
}

bool Directed_Graph::drop_vertex(int u)
{
  if (!Graph::drop_vertex(u)) return false;
  // We need to re-index the poset, which is rather complicated...
  int i,j,nv = nvertex + 1;
  RELATION rho;
  for(i=nv-1; i>=0; --i) {
    for(j=nv-1; j>=0; --j) {
      if (i < u && j < u) continue;
      if (i == j) continue;
      rho = orientation.get_order(i,j);
      if (rho == INCOMPARABLE) continue;
      if (i >= u && j >= u) {
        if (rho == BEFORE) {
          orientation.set_order(i-1,j-1);
        }
        else {
          orientation.set_order(j-1,i-1);
        }
      }
      else if (i >= u) {
        if (rho == BEFORE) {
          orientation.set_order(i-1,j);
        }
        else {
          orientation.set_order(j,i-1);
        }
      }
      else {
        // Has to be j >= u
        if (rho == BEFORE) {
          orientation.set_order(i,j-1);
        }
        else {
          orientation.set_order(j-1,i);
        }
      }
    }
  }
  orientation.N -= 1;
  return true; 
}

bool Directed_Graph::drop_edge(int u,int v)
{
  if (!Graph::drop_edge(u,v)) return false;
  orientation.unset_order(u,v);
  return true;  
}

bool Directed_Graph::foliation_x(int u,int v)
{
  if (u == v) return false;
  if (neighbours[u].empty()) return false;
  if (v == -1) v = RND.irandom(neighbours[u]);
  return drop_edge(u,v);
}

bool Directed_Graph::foliation_m(int u,int v)
{
  if (v == -1) {
    do {
      v = RND.irandom(nvertex);
      if (v != u) break;
    } while(true);
  }
  return add_edge(u,v);
}

int Directed_Graph::fission_x(int v)
{
  int u = nvertex;
  nvertex++;
  add_edge(v,u);
  return u;
}

int Directed_Graph::fission_m(int v)
{
  int u = nvertex;
  std::set<int> S = neighbours[v];
  std::set<int>::const_iterator it;
  nvertex++;
  add_edge(v,u);
  for(it=S.begin(); it!=S.end(); ++it) {
    add_edge(u,*it);
  }
  return u;
}

bool Directed_Graph::add_edge(int u,int v)
{
  if (!Graph::add_edge(u,v)) return false;
  double alpha = RND.drandom();
  if (alpha < 0.1) {
    orientation.set_order(u,v);
  }
  else if (alpha < 0.2) {
    orientation.set_order(v,u);
  }
  return true;
}

bool Directed_Graph::add_edge(int u,int v,RELATION rho)
{
  if (!Graph::add_edge(u,v)) return false;
  if (rho == BEFORE) {
    return orientation.set_order(u,v);
  }
  else if (rho == AFTER) {
    return orientation.set_order(v,u);
  }
  return true;
}

bool Directed_Graph::path_connected(int u,int v) const
{
  // Is it possible to get from u to v following the orientation of the graph edges?
  int i,j;
  bool output = false;
  std::set<int> current,next,S;
  std::set<int>::const_iterator it,jt;
  hash_map::const_iterator qt;

  current.insert(u);
  do {
    for(it=current.begin(); it!=current.end(); ++it) {
      i = *it;
      for(jt=neighbours[i].begin(); jt!=neighbours[i].end(); ++jt) {
        j = *jt;
        S.clear();
        S.insert(i); S.insert(j);
        qt = index_table.find(S);
        if (edges[qt->second].active && orientation.get_order(i,j) == BEFORE) next.insert(j);
      }
    }
    if (next.empty()) break;
    if (next.count(v) > 0) {
      output = true;
      break;
    }
    current = next;
    next.clear();
  } while(true);
  return output;
}

void Directed_Graph::compute_sinks(std::set<int>& output) const
{
  int i,j;
  bool sink;
  std::set<int> S;
  std::set<int>::const_iterator it;
  hash_map::const_iterator qt;

  output.clear();
  // A sink is a vertex all of whose edges are incoming...
  for(i=0; i<nvertex; ++i) {
    if (neighbours[i].empty()) continue;
    sink = true;
    for(it=neighbours[i].begin(); it!=neighbours[i].end(); ++it) {
      j = *it;
      S.clear();
      S.insert(i); S.insert(j);
      qt = index_table.find(S);
      if (edges[qt->second].active && orientation.get_order(i,j) == BEFORE) {
        // If this vertex has an outgoing edge, it isn't a sink...
        sink = false;
        break;
      }
    }
    if (sink) output.insert(i);
  }
}

void Directed_Graph::compute_sources(std::set<int>& output) const
{
  int i,j;
  bool source;
  std::set<int> S;
  std::set<int>::const_iterator it;
  hash_map::const_iterator qt;

  output.clear();
  // A source is a vertex all of whose edges are outgoing...
  for(i=0; i<nvertex; ++i) {
    if (neighbours[i].empty()) continue;
    source = true;
    for(it=neighbours[i].begin(); it!=neighbours[i].end(); ++it) {
      j = *it;
      S.clear();
      S.insert(i); S.insert(j);
      qt = index_table.find(S);
      if (edges[qt->second].active && orientation.get_order(i,j) == AFTER) {
        // If this vertex has an incoming edge, it isn't a source...
        source = false;
        break;
      }
    }
    if (source) output.insert(i);
  }
}
