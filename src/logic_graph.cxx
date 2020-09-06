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
  if (propositional_density < std::numeric_limits<double>::epsilon()) throw std::invalid_argument("The propositional density must be greater than zero!");

  std::set<int> atoms;
  // Make sure this set has at least one member...
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
  if (nvertex > 0) logic = new Propositional_System(*source.logic);
}

Logic_Graph& Logic_Graph::operator =(const Logic_Graph& source)
{
  if (this == &source) return *this;

  clear();

  nvertex = source.nvertex;
  neighbours = source.neighbours;
  edges = source.edges;
  logical_breadth = source.logical_breadth;
  if (nvertex > 0) logic = new Propositional_System(*source.logic);

  return *this;
}

Logic_Graph::~Logic_Graph()
{
  if (nvertex > 0) delete logic;
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

unsigned int Logic_Graph::compute_logical_breadth()
{
  int i,j;
  unsigned int n,ntotal = 0;
  std::set<int>::const_iterator it,jt;
  std::set<int> current,atoms;

  logical_breadth.clear();

  for(i=0; i<nvertex; ++i) {
    // First count the number of atomic propositions in my own
    // proposition...
    logic->get_atoms(i,atoms);
    // Then get the atoms contained in all the propositions of my neighbours...
    for(it=neighbours[i].begin(); it!=neighbours[i].end(); ++it) {
      j = *it;
      logic->get_atoms(j,current);
      for(jt=current.begin(); jt!=current.end(); ++it) {
        atoms.insert(*it);
      }
    }
    n = atoms.size();
    logical_breadth.push_back(n);
    ntotal += n;
  }
  return ntotal;
}

double Logic_Graph::rationalize_topology(const std::string& type)
{
  // We begin this method with a complete graph topology, which is
  // obviously connected. At each iteration, we examine an edge to
  // verify that the propositions of its two vertices are consistent
  // under the AND operator, otherwise we remove the edge (unless 
  // this would disconnect the graph).
  int i,j,vx[2];
  unsigned int m,n,rsum = 0,its = 0;
  std::vector<unsigned int> bcount;
  std::set<int>::const_iterator it;
  const unsigned int N = order()/2;

  make_complete();

  for(i=0; i<nvertex; ++i) {
    bcount.push_back(logic->bit_count(i));
  }

  do {
    m = RND.irandom(edges.size());
    its++;
    edges[m].get_vertices(vx);
    n = logic->consistency(vx[0],vx[1],type);
    // Keep this edge if the bit count is complete...
    if (n == bcount[vx[0]]) continue;
    // If not, drop it...
    drop_edge(vx[0],vx[1]);
    // Unless this would disconnect the graph...
    if (!connected()) add_edge(vx[1],vx[1]);
  } while(its < N);

  for(i=0; i<nvertex; ++i) {
    for(it=neighbours[i].begin(); it!=neighbours[i].end(); ++it) {
      j = *it;
      n = logic->consistency(i,j,type);
      rsum += bcount[i] - n;
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

void Logic_Graph::add_theorem(int v)
{
  int i;
  std::set<int> atoms,current;
  std::set<int>::const_iterator it,jt;

  logic->get_atoms(v,atoms);
  // Then get the atoms contained in all the propositions of its neighbours...
  for(it=neighbours[v].begin(); it!=neighbours[v].end(); ++it) {
    i = *it;
    logic->get_atoms(i,current);
    for(jt=current.begin(); jt!=current.end(); ++it) {
      atoms.insert(*it);
    }
  }
  logic->add_theorem(atoms);
}

int Logic_Graph::fission_x(int v)
{
  unsigned int n = Graph::fission_x(v);

  add_theorem(v);

  return n;
}

int Logic_Graph::fission_m(int v)
{
  unsigned int n = Graph::fission_m(v);

  add_theorem(v);

  return n;
}

