#include "logic_graph.h"

using namespace SYNARMOSMA;

extern Random RND;

Logic_Graph::Logic_Graph() : Graph()
{

}

Logic_Graph::Logic_Graph(int n,double propositional_density) : Graph(n,true)
{
  // We will choose a number of atoms that is a multiple 
  // of the number of vertices...
  assert(propositional_density > std::numeric_limits<double>::epsilon());

  std::set<int> atoms;
  // Make sure the set isn't empty...
  atoms.insert(0);
  for(int i=1; i<nvertex; ++i) {
    if (RND.drandom() < propositional_density) atoms.insert(i);
  }
  logic = new Propositional_System(atoms,nvertex);
}

Logic_Graph::Logic_Graph(const Logic_Graph& source)
{
  nvertex = source.nvertex;
  neighbours = source.neighbours;
  edges = source.edges;
  logical_breadth = source.logical_breadth;
  if (nvertex > 0) {
    logic = new Propositional_System;
    logic = source.logic;
  }
}

Logic_Graph& Logic_Graph::operator =(const Logic_Graph& source)
{
  if (this == &source) return *this;

  clear();

  nvertex = source.nvertex;
  neighbours = source.neighbours;
  edges = source.edges;
  logical_breadth = source.logical_breadth;
  if (nvertex > 0) {
    logic = new Propositional_System;
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
  if (nvertex > 0) delete logic;
  nvertex = 0;
  edges.clear();
  index_table.clear();
  neighbours.clear();
}

int Logic_Graph::serialize(std::ofstream& s) const
{
  int i,j,count = 0,nedge = (signed) edges.size();
  std::set<int>::const_iterator it;

  s.write((char*)(&nvertex),sizeof(int)); count += sizeof(int);
  for(i=0; i<nvertex; ++i) {
    j = (signed) neighbours[i].size();
    s.write((char*)(&j),sizeof(int)); count += sizeof(int);
    for(it=neighbours[i].begin(); it!=neighbours[i].end(); ++it) {
      j = *it;
      s.write((char*)(&j),sizeof(int)); count += sizeof(int);
    }
  }
  s.write((char*)(&nedge),sizeof(int)); count += sizeof(int);
  for(i=0; i<nedge; ++i) {
    count += edges[i].serialize(s);
  }
  count += logic->serialize(s);

  return count;
}

int Logic_Graph::deserialize(std::ifstream& s)
{
  int i,j,k,n,vx[2],count = 0;
  Edge q;
  std::set<int> S;

  clear();

  s.read((char*)(&nvertex),sizeof(int)); count += sizeof(int);
  for(i=0; i<nvertex; ++i) {
    s.read((char*)(&n),sizeof(int)); count += sizeof(int);
    for(j=0; j<n; ++j) {
      s.read((char*)(&k),sizeof(int)); count += sizeof(int);
      S.insert(k);
    }
    neighbours.push_back(S);
    S.clear();
  }
  s.read((char*)(&n),sizeof(int)); count += sizeof(int);
  for(i=0; i<n; ++i) {
    count += q.deserialize(s);
    edges.push_back(q);
    q.get_vertices(vx);
    S.clear();
    S.insert(vx[0]); S.insert(vx[1]);
    index_table[S] = i;
  }
  count += logic->deserialize(s);
  compute_logical_breadth();

  return count;
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
    logic->get_atoms(i,atoms);
    for(it=neighbours[i].begin(); it!=neighbours[i].end(); ++it) {
      in1 = *it;
      logic->get_atoms(in1,current);
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
  std::string op = "AND";
  std::vector<int> bcount;
  std::set<int>::const_iterator it;
  const int N = order()/2;

  make_complete();

  for(i=0; i<nvertex; ++i) {
    bcount.push_back(logic->bit_count(i));
  }

  do {
    m = RND.irandom(edges.size());
    its++;
    edges[m].get_vertices(vx);
    n = logic->consistency(vx[0],vx[1],op);
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
      n = logic->consistency(i,m,op);
      rsum += (bcount[i] - n);
    }
  }
  return double(rsum)/double(nvertex);
}

bool Logic_Graph::drop_vertex(int v)
{
  if (Graph::drop_vertex(v)) {
    if (logic->drop_theorem(v)) return true;
  }
  return false;
}

bool Logic_Graph::fusion(int v,int u)
{
  if (Graph::fusion(v,u)) {
    if (logic->drop_theorem(u)) return true;
  }
  return false;
}

int Logic_Graph::fission_x(int v)
{
  unsigned int i,n = Graph::fission_x(v),na = logic->get_number_atoms();
  std::set<int> atoms;

  for(i=0; i<na; ++i) {
    atoms.insert(i);
  }
  logic->add_theorem(atoms);
  return n;
}

int Logic_Graph::fission_m(int v)
{
  unsigned int i,n = Graph::fission_m(v),na = logic->get_number_atoms();
  std::set<int> atoms;

  for(i=0; i<na; ++i) {
    atoms.insert(i);
  }
  logic->add_theorem(atoms);
  return n;
}

