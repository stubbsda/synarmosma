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
  int i,j,nc = 0;
  double alpha;
  std::set<int> vx;
  RELATION rho;

  for(i=0; i<n; ++i) {
    for(j=1+i; j<n; ++j) {
      vx.insert(i);
      vx.insert(j);
      index_table[vx] = nc;
      rho = DISPARATE;
      alpha = RND.drandom();
      if (alpha < 0.15) {
        rho = BEFORE;
      }
      else if (alpha < 0.3) {
        rho = AFTER;
      }      
      edges.push_back(Edge(i,j,rho));
      neighbours[i].insert(j);
      neighbours[j].insert(i);
      vx.clear();
      nc++;
    }
  }
}

Directed_Graph::Directed_Graph(int n,double p) : Graph(n)
{
  int i,j,nc = 0;
  double alpha;
  std::set<int> vx;
  RELATION rho;

  for(i=0; i<n; ++i) {
    for(j=1+i; j<n; ++j) {
      alpha = RND.drandom();
      if (alpha > p) continue;
      vx.insert(i);
      vx.insert(j);
      index_table[vx] = nc;
      rho = DISPARATE;
      alpha = RND.drandom();
      if (alpha < 0.15) {
        rho = BEFORE;
      }
      else if (alpha < 0.3) {
        rho = AFTER;
      }      
      edges.push_back(Edge(i,j,rho));
      neighbours[i].insert(j);
      neighbours[j].insert(i);
      vx.clear();
      nc++;
    }
  }
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

int Directed_Graph::distance(int u,int v) const
{
  // A method to calculate the topological distance from 
  // u to v; the method returns -1 if it is impossible to 
  // get from u to v.
  if (u == v) return 0;

  int i,j,l = -1,its = 1;
  bool visited[nvertex];
  RELATION rho;
  std::set<int> S,current,next;
  std::set<int>::const_iterator it,jt;
  hash_map::const_iterator qt;

  for(i=0; i<nvertex; ++i) {
    visited[i] = false;
  }
  current.insert(u);
  visited[u] = true;
  do {
    for(it=current.begin(); it!=current.end(); ++it) {
      i = *it;
      for(jt=neighbours[i].begin(); jt!=neighbours[i].end(); ++jt) {
        j = *jt;
        if (visited[j]) continue;
        S.clear();
        S.insert(i); S.insert(j);
        qt = index_table.find(S);
        if (!edges[qt->second].active) continue;
        rho = edges[qt->second].direction;
        if (i < j && rho == BEFORE) {
          next.insert(j);
        }
        else if (i > j && rho == AFTER) {
          next.insert(j);
        }
      }
    }
    if (next.empty()) break;
    if (next.count(v) > 0) {
      l = its;
      break;
    }
    for(it=next.begin(); it!=next.end(); ++it) {
      visited[*it] = true;
    }
    current = next;
    next.clear();
    its++;
  } while(true);
  return l;
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
  if (rho == DISPARATE) { 
    edges[qt->second].direction = rho;
  }
  else {
    if (u < v) {
      edges[qt->second].direction = rho;
    }
    else {
      edges[qt->second].direction = (rho == BEFORE) ? AFTER : BEFORE;
    }
  }
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
  RELATION rho;
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
        if (!edges[qt->second].active) continue;
        rho = edges[qt->second].direction;
        if ((i < j && rho == BEFORE) || (j < i && rho == AFTER)) next.insert(j);
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

bool Directed_Graph::directed_cycle(const std::vector<int>& path,int base,int length) const
{
  if (base == path[0] && path.size() > 2) return true;
  if ((signed) path.size() == length) return false;
  int v;
  bool out;
  std::set<int> S;
  std::set<int>::const_iterator it;
  RELATION rho;
  hash_map::const_iterator qt;

  for(it=neighbours[base].begin(); it!=neighbours[base].end(); ++it) {
    v = *it;
    S.clear();
    S.insert(base); S.insert(v);
    qt = index_table.find(S);
    if (!edges[qt->second].active) continue;
    rho = edges[qt->second].get_direction();
    if ((base < v && rho == BEFORE) || (v < base && rho == AFTER)) {
      std::vector<int> npath = path;
      npath.push_back(v);
      out = directed_cycle(npath,v,length);
      if (out) return true;
    }
  }
  return false;
}

bool Directed_Graph::acyclic() const
{
  std::vector<int> path;

  for(int i=0; i<nvertex; ++i) {
    path.clear();
    path.push_back(i);
    if (directed_cycle(path,i,nvertex)) return false;
  }
  return true;
}

void Directed_Graph::compute_sinks(std::set<int>& output) const
{
  int i,j;
  bool sink;
  std::set<int> S;
  std::set<int>::const_iterator it;
  RELATION rho;
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
      if (!edges[qt->second].active) continue;
      rho = edges[qt->second].direction;
      if ((i < j && rho == BEFORE) || (j < i && rho == AFTER)) {
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
  RELATION rho;
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
      if (!edges[qt->second].active) continue;
      rho = edges[qt->second].direction;
      if ((i < j && rho == AFTER) || (j < i && rho == BEFORE)) {
        // If this vertex has an incoming edge, it isn't a source...
        source = false;
        break;
      }
    }
    if (source) output.insert(i);
  }
}
