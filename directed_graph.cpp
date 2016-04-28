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
    if (it->direction == DISPARATE) null++;
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

void Directed_Graph::compute_distances(edge_hash& output) const
{
  int i,j,k,delta;
  RELATION rho;
  std::pair<int,int> pr;
  std::vector<int> distances;
  std::set<int> S;
  std::set<int>::const_iterator it;
  hash_map::const_iterator qt;

  for(i=0; i<nvertex; ++i) {
    for(j=0; j<nvertex; ++j) {
      distances.push_back(std::numeric_limits<int>::max()); 
    }
  }
  for(i=0; i<nvertex; ++i) {
    distances[nvertex*i+i] = 0;
    for(it=neighbours[i].begin(); it!=neighbours[i].end(); ++it) {
      j = *it;
      S.clear();
      S.insert(i); S.insert(j);
      qt = index_table.find(S);
      rho = edges[qt->second].direction;
      if (i < j && rho == BEFORE) {      
        distances[nvertex*i+j] = 1;
      }
      else if (i > j && rho == AFTER) {
        distances[nvertex*i+j] = 1;
      }
    }
  }

  for(k=0; k<nvertex; ++k) {
    for(i=0; i<nvertex; ++i) {
      if (distances[nvertex*i+k] == std::numeric_limits<int>::max()) continue;
      for(j=0; j<nvertex; ++j) {
        if (distances[nvertex*k+j] == std::numeric_limits<int>::max()) continue;
        delta = distances[nvertex*i+k] + distances[nvertex*k+j];
        if (delta < distances[nvertex*i+j]) distances[nvertex*i+j] = delta;
      }
    }
  }
  output.clear();
  for(i=0; i<nvertex; ++i) {
    for(j=0; j<nvertex; ++j) {
      if (i == j) continue;
      if (distances[nvertex*i+j] == std::numeric_limits<int>::max()) distances[nvertex*i+j] = -1;
      pr.first = i; pr.second = j;
      output[pr] = distances[nvertex*i+j];
      pr.first = j; pr.second = i;
      output[pr] = distances[nvertex*j+i];
    }
  }
}

bool Directed_Graph::add_edge(int u,int v,double ell,RELATION rho)
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
  edges[qt->second].capacity = ell;
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
    rho = edges[qt->second].direction;
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

double Directed_Graph::compute_flow(int source,int sink) 
{
  int i;
  bool valid = false;
  RELATION rho;
  std::set<int> S;
  std::set<int>::const_iterator it;
  hash_map::const_iterator qt;
  const int ne = size();

  for(i=0; i<ne; ++i) {
    if (edges[i].direction == DISPARATE) edges[i].capacity = 0.0;
  }
  // Next we need to verify that there is at least one outgoing edge with capacity > 0 
  // for the source and at least one incoming edge with capacity > 0 for the sink
  for(it=neighbours[source].begin(); it!=neighbours[source].end(); ++it) {
    i = *it;
    S.clear();
    S.insert(source);
    S.insert(i);
    qt = index_table.find(S);
    if (edges[qt->second].capacity < std::numeric_limits<double>::epsilon()) continue;
    rho = edges[qt->second].direction;
    if ((source < i && rho == BEFORE) || (source > i && rho == AFTER)) valid = true;
    if (valid) break;
  }
  if (!valid) {
    for(i=0; i<ne; ++i) {
      edges[i].flow = 0.0;
    }
    return 0.0;
  }
  // Now the sink
  valid = false;
  for(it=neighbours[sink].begin(); it!=neighbours[sink].end(); ++it) {
    i = *it;
    S.clear();
    S.insert(sink);
    S.insert(i);
    qt = index_table.find(S);
    if (edges[qt->second].capacity < std::numeric_limits<double>::epsilon()) continue;
    rho = edges[qt->second].direction;
    if ((sink < i && rho == AFTER) || (sink > i && rho == BEFORE)) valid = true;
    if (valid) break;
  }
  if (!valid) {
    for(i=0; i<ne; ++i) {
      edges[i].flow = 0.0;
    }
    return 0.0;
  }
  // So we can finally proceed to computing a non-trivial flow on this directed graph
  int vx[2];
  std::pair<int,int> pr;
  edge_hash rgraph;
  edge_hash::const_iterator qtt;
  // Create a residual graph and fill the residual graph with
  // given capacities in the original graph as residual capacities
  // in residual graph
  for(i=0; i<ne; ++i) {
    edges[i].get_vertices(vx);
    if (edges[i].direction == BEFORE) {
      pr.first = vx[0]; pr.second = vx[1];
      rgraph[pr] = int(10000.0*edges[i].capacity);
    }
    else if (edges[i].direction == AFTER) {
      pr.first = vx[1]; pr.second = vx[0];
      rgraph[pr] = int(10000.0*edges[i].capacity);
    }
  }

  int max_flow = network_flow(rgraph,source,sink,nvertex);
  for(i=0; i<ne; ++i) {
    edges[i].flow = 0.0;
    edges[i].get_vertices(vx);
    if (edges[i].direction == BEFORE) {
      pr.first = vx[0]; pr.second = vx[1]; 
      edges[i].flow = edges[i].capacity - double(rgraph[pr])/10000.0;
    }
    else if (edges[i].direction == AFTER) {
      pr.first = vx[1]; pr.second = vx[0];
      edges[i].flow = edges[i].capacity - double(rgraph[pr])/10000.0;
    }
  }
  return double(max_flow)/10000.0;
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
