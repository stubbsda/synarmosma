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

Logic_Graph::~Logic_Graph()
{
  delete logic;
}

void Logic_Graph::create()
{
  rationalize_topology();
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

int Logic_Graph::rationalize_topology()
{
  // We begin this method with a random graph topology, which is
  // however connected. At each iteration, we examine the causal
  // edges to verify that the proposition of the antecedent vertex
  // implies the propositions of its descendent vertex neighbours
  int i,j,k,tt,test,in1,in2,nfalse,its = 0;
  bool flag,good;
  std::set<int> S;
  std::set<int>::const_iterator it;
  int* ntrue = new int[nvertex];
  std::vector<int>* false_edges = new std::vector<int>[nvertex];
  const int max_iter = 500;

  for(i=0; i<nvertex; ++i) {
    false_edges[i].reserve(100);
  }

  for(i=0; i<nvertex; ++i) {
    ntrue[i] = logic->bit_count(i);
  }

  do {
    its++;
    nedge = 0;
    for(i=0; i<nvertex; ++i) {
      false_edges[i].clear();
    }
    for(i=0; i<nvertex; ++i) {
      for(it=neighbours[i].begin(); it!=neighbours[i].end(); ++it) {
        in1 = *it;
        if (in1 > i) continue;
        nedge += 1;
        in2 = logic->consistency(i,in1,"and");
        if (in2 < ntrue[i]) {
          false_edges[i].push_back(j);
          for(it=neighbours[in1].begin(); it!=neighbours[in1].end(); ++it) {
            if (*it == i) {
              false_edges[in1].push_back(k);
              break;
            }
          }
        }
      }
    }
    nfalse = 0;
    for(i=0; i<nvertex; ++i) {
      nfalse += false_edges[i].size();
    }
    nfalse /= 2;
    if (nfalse == 0 || its == max_iter) break;
    for(i=0; i<nvertex; ++i) {
      if (false_edges[i].empty()) continue;
      tt = false_edges[i].size();
      for(j=0; j<tt; ++j) {
        in1 = false_edges[i][j];
        for(it=neighbours[i].begin(); it!=neighbours[i].end(); ++it) {
          if (*it == in1) continue;
          S.insert(*it);
        }
        neighbours[i] = S;
        S.clear();
      }
      for(j=0; j<tt/2; ++j) {
        flag = true;
        do {
          test = RND.irandom(nvertex);
          if (test == i) continue;
          good = (neighbours[i].count(test) == 1) ? false : true;
          if (!good) continue;
          neighbours[i].insert(test);
          neighbours[test].insert(i);
          flag = false;
        } while(flag);
      }
    }
    if (!connected()) {
      std::cerr << "Logic_Graph instance disconnected in topology rationalization!" << std::endl;
      std::exit(1);
    }
  } while(true);

  // Free the memory
  delete[] false_edges;
  delete[] ntrue;

  return nfalse;
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

