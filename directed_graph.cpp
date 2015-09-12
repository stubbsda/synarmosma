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

Directed_Graph::Directed_Graph(unsigned int n) : Graph(n)
{
  orientation = new Poset(n);
}

Directed_Graph::~Directed_Graph()
{
  delete orientation;
}

Directed_Graph::Directed_Graph(const Directed_Graph& source)
{
  vertices = source.vertices;
  edges = source.edges;
  orientation = new Poset(*source.orientation);
}

Directed_Graph& Directed_Graph::operator =(const Directed_Graph& source) 
{
  if (this == &source) return *this;
  vertices = source.vertices;
  edges = source.edges;
  orientation = new Poset(*source.orientation);
  return *this;
}

void Directed_Graph::clear()
{
  vertices.clear();
  edges.clear();
  orientation->clear();
}

void Directed_Graph::serialize(std::ofstream& s) const
{
  unsigned int i,j,order = get_order();
  bool active;
  std::set<unsigned int>::const_iterator it;

  s.write((char*)(&order),sizeof(int));
  for(i=0; i<order; ++i) {
    active = vertices[i].first;
    s.write((char*)(&active),sizeof(bool));
    j = vertices[i].second.size();
    s.write((char*)(&j),sizeof(int));
    for(it=vertices[i].second.begin(); it!=vertices[i].second.end(); ++it) {
      j = *it;
      s.write((char*)(&j),sizeof(int));
    }
  }
  j = edges.size();
  s.write((char*)(&j),sizeof(int));
  for(i=0; i<j; ++i) {
    edges[i].serialize(s);
  }
  orientation->serialize(s);
}

void Directed_Graph::deserialize(std::ifstream& s)
{
  unsigned int i,j,k,n,m;
  bool active;
  std::set<unsigned int> S;
  Edge q;

  clear();

  s.read((char*)(&n),sizeof(int));
  orientation = new Poset(n);
  for(i=0; i<n; ++i) {
    s.read((char*)(&active),sizeof(bool));
    s.read((char*)(&m),sizeof(int));
    for(j=0; j<m; ++j) {
      s.read((char*)(&k),sizeof(int));
      S.insert(k);
    }
    vertices.push_back(std::pair<bool,std::set<unsigned int> >(active,S));
    S.clear();
  }
  s.read((char*)(&n),sizeof(int));
  for(i=0; i<n; ++i) {
    q.deserialize(s);
    edges.push_back(q);
  }
  orientation->deserialize(s);
}

bool Directed_Graph::foliation_x(unsigned int v)
{
  if (!vertices[v].first) return false;
  unsigned int w,order = get_order();
  do {
    w = RND.irandom(order);
    if (vertices[w].first && w != v) break;
  } while(true);
  return foliation_x(v,w);
}

bool Directed_Graph::foliation_x(unsigned int v1,unsigned int v2)
{
  if (v1 == v2) return false;
  if (!vertices[v1].first || vertices[v1].second.empty()) return false;
  std::set<unsigned int> S;

  S.insert(v1); S.insert(v2);
  vertices[v1].second.erase(v2);
  vertices[v2].second.erase(v1);
  hash_map::const_iterator qt = index_table.find(S);
  edges[qt->second].active = false;
  return true;
}

bool Directed_Graph::foliation_m(unsigned int v)
{
  if (!vertices[v].first) return false;
  unsigned int w,order = get_order();
  do {
    w = RND.irandom(order);
    if (vertices[w].first && w != v) break;
  } while(true);
  return foliation_m(v,w);
}

bool Directed_Graph::foliation_m(unsigned int v1,unsigned int v2)
{
  if (!vertices[v1].first || !vertices[v1].first) return false;
  return add_edge(v1,v2);
}

unsigned int Directed_Graph::fission_x(unsigned int v)
{
  assert(vertices[v].first);
  unsigned int u = add_vertex();
  add_edge(v,u);
  return u;
}

unsigned int Directed_Graph::fission_m(unsigned int v)
{
  unsigned int u = add_vertex();
  std::set<unsigned int> S = vertices[v].second;
  std::set<unsigned int>::const_iterator it;

  add_edge(v,u);
  for(it=S.begin(); it!=S.end(); ++it) {
    add_edge(u,*it);
  }
  return u;
}

unsigned int Directed_Graph::add_vertex()
{  
  orientation->add_element();
  return Graph::add_vertex();
}

bool Directed_Graph::drop_vertex(unsigned int v)
{
  if (!Graph::drop_vertex(v)) return false;
  orientation->drop_element(v);
  return true;
}

bool Directed_Graph::add_edge(unsigned int v1,unsigned int v2)
{
  if (v1 == v2) return false;
  if (!Graph::add_edge(v1,v2)) return false;

  // Choose a random orientation for this edge...
  double alpha = RND.drandom();
  if (alpha < 0.25) { 
    orientation->set_order(v1,v2);
  }
  else if (alpha < 0.5) {

  }
  else {

  }
  return true;
}

bool Directed_Graph::add_edge(unsigned int v1,unsigned int v2,RELATION rho)
{
  if (v1 == v2) return false;
  if (!Graph::add_edge(v1,v2)) return false;

  if (v1 < v2) {

  }
  else {

  }

  return true;
}

bool Directed_Graph::drop_edge(unsigned int v1,unsigned int v2)
{
  if (!Graph::drop_edge(v1,v2)) return false;

  orientation->unset_order(v1,v2);

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
