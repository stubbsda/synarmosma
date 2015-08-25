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
  nedge = source.nedge;
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
  nedge = source.nedge;
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

void Logic_Graph::serialize(std::ofstream& s) const
{

}

void Logic_Graph::deserialize(std::ifstream& s)
{

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
  int i,p,q,n,rsum,its = 0;
  std::vector<int> bcount;
  std::set<int>::const_iterator it;

  make_complete();

  for(i=0; i<nvertex; ++i) {
    bcount.push_back(logic->bit_count(i));
  }

  do {
    p = RND.irandom(nvertex);
    q = RND.irandom(nvertex);
    if (p == q) continue;
    its++;
    if (!connected(p,q)) continue;
    n = logic->consistency(p,q,"and");
    // Keep this edge if the bit count is complete...
    if (n == bcount[p]) continue;
    // If not, drop it
    drop_edge(p,q);
    // Unless this would disconnect the graph...
    if (!connected()) add_edge(p,q);
  } while(its < nvertex*nvertex);

  rsum = 0;
  for(i=0; i<nvertex; ++i) {
    for(it=neighbours[i].begin(); it!=neighbours[i].end(); ++it) {
      p = *it;
      n = logic->consistency(i,p,"and");
      rsum += (bcount[i] - n);
    }
  }
  return double(rsum)/double(nvertex);
}

bool Logic_Graph::amputation(int v)
{
  if (Graph::amputation(v)) {
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

