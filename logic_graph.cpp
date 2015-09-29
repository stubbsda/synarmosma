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

#include "logic_graph.h"

using namespace SYNARMOSMA;

extern Random RND;

Logic_Graph::Logic_Graph() : Graph()
{

}

Logic_Graph::Logic_Graph(int n) : Graph(n)
{
  logic = new Propositional_System(nvertex);
}

Logic_Graph::Logic_Graph(const Logic_Graph& source)
{
  nvertex = source.nvertex;
  neighbours = source.neighbours;
  edges = source.edges;
  logical_breadth = source.logical_breadth;
  if (nvertex > 0) {
    logic = new Propositional_System(nvertex);
    logic = source.logic;
  }
}

Logic_Graph& Logic_Graph::operator =(const Logic_Graph& source)
{
  if (this == &source) return *this;
  if (nvertex > 0) delete logic;
  nvertex = source.nvertex;
  neighbours = source.neighbours;
  edges = source.edges;
  logical_breadth = source.logical_breadth;
  if (nvertex > 0) {
    logic = new Propositional_System(nvertex);
    logic = source.logic;
  }
  return *this;
}

Logic_Graph::~Logic_Graph()
{
  if (nvertex > 0) delete logic;
}

void Logic_Graph::create()
{
  rationalize_topology();
}

void Logic_Graph::clear()
{
  logic->clear();
  nvertex = 0;
  edges.clear();
  index_table.clear();
  neighbours.clear();
}

void Logic_Graph::serialize(std::ofstream& s) const
{
  int i,j,nedge = edges.size();
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
  s.write((char*)(&nedge),sizeof(int));
  for(i=0; i<nedge; ++i) {
    edges[i].serialize(s);
  }
  logic->serialize(s);
}

void Logic_Graph::deserialize(std::ifstream& s)
{
  int i,j,k,n;
  Edge q;
  std::set<int> S;

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
    index_table[q.nodes] = i;
  }
  logic->deserialize(s);
  compute_logical_breadth();
}

void Logic_Graph::compute_logical_breadth()
{
  int i,in1;
  std::set<int>::const_iterator it,jt;
  std::set<int> current,atoms;

  logical_breadth.clear();

  for(i=0; i<nvertex; ++i) {
    // First count the number of atomic propositions in my own
    // proposition:
    logic->theorems[i].atoms(atoms);
    for(it=neighbours[i].begin(); it!=neighbours[i].end(); ++it) {
      in1 = *it;
      logic->theorems[in1].atoms(current);
      for(jt=current.begin(); jt!=current.end(); ++it) {
        atoms.insert(*it);
      }
    }
    logical_breadth.push_back(atoms.size());
  }
}

double Logic_Graph::rationalize_topology()
{
  // We begin this method with a complete graph topology, which is
  // obviously connected. At each iteration, we examine the causal
  // edges to verify that the proposition of the antecedent vertex
  // implies the propositions of its descendent vertex neighbours
  int i,m,n,vx[2],rsum,its = 0;
  std::vector<int> bcount;
  std::set<int>::const_iterator it;
  const int N = order()/2;

  make_complete();

  for(i=0; i<nvertex; ++i) {
    bcount.push_back(logic->bit_count(i));
  }

  do {
    m = RND.irandom(edges.size());
    if (!edges[m].active) continue;
    its++;
    edges[m].get_nodes(vx);
    n = logic->consistency(vx[0],vx[1],"and");
    // Keep this edge if the bit count is complete...
    if (n == bcount[vx[0]]) continue;
    // If not, drop it
    drop_edge(vx[0],vx[1]);
    // Unless this would disconnect the graph...
    if (!connected()) add_edge(vx[1],vx[1]);
  } while(its < N);

  rsum = 0;
  for(i=0; i<nvertex; ++i) {
    for(it=neighbours[i].begin(); it!=neighbours[i].end(); ++it) {
      m = *it;
      n = logic->consistency(i,m,"and");
      rsum += (bcount[i] - n);
    }
  }
  return double(rsum)/double(nvertex);
}

bool Logic_Graph::drop_vertex(int v)
{
  if (Graph::drop_vertex(v)) {
    logic->theorems.erase(logic->theorems.begin() + v);
    logic->truth.erase(logic->truth.begin() + v);
    return true;
  }
  return false;
}

bool Logic_Graph::foliation_m(int v,int u)
{
  return Graph::foliation_m(v,u);
}

bool Logic_Graph::foliation_x(int v,int u)
{
  return Graph::foliation_x(v,u);
}

bool Logic_Graph::add_edge(int v,int u)
{
  return Graph::add_edge(v,u);
}

bool Logic_Graph::fusion(int v,int u)
{
  if (Graph::fusion(v,u)) {
    logic->theorems.erase(logic->theorems.begin() + u);
    logic->truth.erase(logic->truth.begin() + u);
    return true;
  }
  return false;
}

int Logic_Graph::fission_x(int v)
{
  int n = Graph::fission_x(v);
  logic->theorems.push_back(Proposition(logic->natom));
  logic->truth.push_back(boost::dynamic_bitset<>(logic->nuniverse));
  return n;
}

int Logic_Graph::fission_m(int v)
{
  int n = Graph::fission_m(v);
  logic->theorems.push_back(Proposition(logic->natom));
  logic->truth.push_back(boost::dynamic_bitset<>(logic->nuniverse));
  return n;
}

