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
}

Directed_Graph& Directed_Graph::operator =(const Directed_Graph& source) 
{
  if (this == &source) return *this;
  nvertex = source.nvertex;
  neighbours = source.neighbours;
  edges = source.edges;
  index_table = source.index_table;
  return *this;
}

int Directed_Graph::directedness() const
{
  // A method to calculate how many of the edges are directed...
  int output,null = 0;
  std::vector<Edge>::const_iterator it;

  for(it=edges.begin(); it!=edges.end(); ++it) {
    if (it->active) {
      if (it->direction == DISPARATE) null++;
    }
  }
  output = size() - null;
  return output;
}

bool Directed_Graph::add_edge(int u,int v)
{
  if (!Graph::add_edge(u,v)) return false;
  RELATION rho = DISPARATE;
  double alpha = RND.drandom();
  std::set<int> S;
  S.insert(u); S.insert(v);
  hash_map::const_iterator qt = index_table.find(S);
  if (alpha < 0.15) {
    rho = BEFORE;
  }
  else if (alpha < 0.3) {
    rho = AFTER;
  }
  edges[qt->second].direction = rho;
  return true;
}

bool Directed_Graph::add_edge(int u,int v,RELATION rho)
{
  if (!Graph::add_edge(u,v)) return false;
  std::set<int> S;
  S.insert(u); S.insert(v);
  hash_map::const_iterator qt = index_table.find(S);
  edges[qt->second].direction = rho;
  return true;
}

bool Directed_Graph::mutate_edge(int u,int v)
{
  if (u == v) return false;
  std::set<int> S;
  S.insert(u); S.insert(v);
  hash_map::const_iterator qt = index_table.find(S);
  if (qt == index_table.end()) return false;
  int n = qt->second;
  if (!edges[n].active) return false;
  double alpha = RND.drandom();
  RELATION rho = edges[n].direction;
  if (rho == DISPARATE) {
    edges[n].direction = (alpha < 0.5) ? BEFORE : AFTER;
  }
  else {
    if (rho == BEFORE) {
      edges[n].direction = (alpha < 0.25) ? AFTER : DISPARATE;
    }
    else {
      edges[n].direction = (alpha < 0.25) ? BEFORE : DISPARATE;
    }
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
        if (edges[qt->second].active && edges[qt->second].direction == BEFORE) next.insert(j);
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
      if (edges[qt->second].active && edges[qt->second].direction == BEFORE) {
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
      if (edges[qt->second].active && edges[qt->second].direction == AFTER) {
        // If this vertex has an incoming edge, it isn't a source...
        source = false;
        break;
      }
    }
    if (source) output.insert(i);
  }
}
