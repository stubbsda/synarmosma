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

Directed_Graph::Directed_Graph() : Schema() 
{

}

Directed_Graph::Directed_Graph(int n) : Schema(n)
{

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

void Directed_Graph::clear()
{
  nvertex = 0;
  neighbours.clear();
  edges.clear();
  index_table.clear();
}

bool Directed_Graph::foliation_x(int v1,int v2)
{
  if (v1 == v2) return false;
  if (neighbours[v1].empty()) return false;
  if (v2 == -1) v2 = RND.irandom(neighbours[v1]);
  neighbours[v1].erase(v2);
  neighbours[v2].erase(v1);
  std::string clef = make_key(v1,v2);
  string_hash::const_iterator qt = index_table.find(clef + ":-1");
  if (qt != index_table.end()) {
    edges[qt->second].active = false;
    return true;
  }
  qt = index_table.find(clef+":+1");
  edges[qt->second].active = false;
  return true;  
}

bool Directed_Graph::foliation_m(int v1,int v2)
{
  if (v2 == -1) {
    do {
      v2 = RND.irandom(nvertex);
      if (v2 != v1) break;
    } while(true);
  }
  return add_edge(v1,v2);
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

bool Directed_Graph::add_edge(int v1,int v2)
{
  if (v1 == v2) return false;
  bool frwd,bkwd;
  std::string name,base = make_key(v1,v2);
  string_hash::const_iterator qt;

  name = base + ":1";
  qt = index_table.find(name);
  frwd = (qt != index_table.end()) ? true : false;

  name = base + ":-1";
  qt = index_table.find(name);
  bkwd = (qt != index_table.end()) ? true : false;  

  if (frwd && bkwd) return false;

  neighbours[v1].insert(v2);
  neighbours[v2].insert(v1);

  if (frwd) {
    // Add the backward-oriented edge
    edges.push_back(Edge(v1,v2,BACKWARD));
    index_table[base+":-1"] = (signed) edges.size() - 1;
  }
  else if (bkwd) {
    // Add the forward-oriented edge
    edges.push_back(Edge(v1,v2,FORWARD));
    index_table[base+":1"] = (signed) edges.size() - 1;
  }
  else {
    // Randomly choose a direction...
    if (RND.irandom(2) == 0) {
      edges.push_back(Edge(v1,v2,FORWARD));
      index_table[base+":1"] = (signed) edges.size() - 1;
    }
    else {
      edges.push_back(Edge(v1,v2,BACKWARD));
      index_table[base+":-1"] = (signed) edges.size() - 1;
    }
  }
  return true;
}

bool Directed_Graph::add_edge(int v1,int v2,DIRECTION d)
{
  if (v1 == v2) return false;
  string_hash::const_iterator qt;
  std::string name = make_key(v1,v2) + ":";
  if (v1 < v2) {
    name += (d == FORWARD) ? "1" : "-1";
  }
  else {
    name += (d == FORWARD) ? "-1" : "1";
  }
  qt = index_table.find(name);
  if (qt != index_table.end()) return false;
  neighbours[v1].insert(v2);
  neighbours[v2].insert(v1);
  edges.push_back(Edge(v1,v2,d));
  index_table[name] = (signed) edges.size() - 1;
  return true;
}

bool Directed_Graph::path_connected(int u,int v) const
{
  // Is it possible to get from u to v following the orientation of the graph edges?
  int i,j;
  bool output = false;
  std::set<int> current,next;
  std::set<int>::const_iterator it,jt;
  string_hash::const_iterator qt;

  current.insert(u);
  do {
    for(it=current.begin(); it!=current.end(); ++it) {
      i = *it;
      for(jt=neighbours[i].begin(); jt!=neighbours[i].end(); ++jt) {
        j = *jt;
        if (j < i) continue;
        qt = index_table.find(make_key(i,j) + ":1");
        if (qt == index_table.end()) continue;
        if (edges[qt->second].active) next.insert(j);
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

int Directed_Graph::two_cycles() const
{
  // Return how many two-cycles there are in this directed graph
  int i,j,t_cycle = 0;
  std::string key;
  std::set<int>::const_iterator it;
  string_hash::const_iterator qt;

  for(i=0; i<nvertex; ++i) {
    for(it=neighbours[i].begin(); it!=neighbours[i].end(); ++it) {
      j = *it;
      if (j < i) continue;
      key = make_key(i,j);
      qt = index_table.find(key + ":1");
      if (qt == index_table.end()) continue;
      if (!edges[qt->second].active) continue;
      // Now see if there is an edge going back from j to i...
      qt = index_table.find(key + ":-1");
      if (qt == index_table.end()) continue;
      if (!edges[qt->second].active) continue;
      t_cycle++;
    }
  }
  return t_cycle;
}

void Directed_Graph::compute_sinks(std::set<int>& output) const
{
  int i,j;
  bool sink;
  std::string key;
  std::set<int>::const_iterator it;
  string_hash::const_iterator qt;

  output.clear();
  // A sink is a vertex all of whose edges are incoming...
  for(i=0; i<nvertex; ++i) {
    if (neighbours[i].empty()) continue;
    sink = true;
    for(it=neighbours[i].begin(); it!=neighbours[i].end(); ++it) {
      j = *it;
      if (j < i) continue;
      key = make_key(i,j);
      qt = index_table.find(key + ":1");
      if (qt != index_table.end()) {
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
  std::string key;
  std::set<int>::const_iterator it;
  string_hash::const_iterator qt;

  output.clear();
  // A source is a vertex all of whose edges are outgoing...
  for(i=0; i<nvertex; ++i) {
    if (neighbours[i].empty()) continue;
    source = true;
    for(it=neighbours[i].begin(); it!=neighbours[i].end(); ++it) {
      j = *it;
      if (j < i) continue;
      key = make_key(i,j);
      qt = index_table.find(key + ":-1");
      if (qt != index_table.end()) {
        // If this vertex has an incoming edge, it isn't a source...
        source = false;
        break;
      }
    }
    if (source) output.insert(i);
  }
}
