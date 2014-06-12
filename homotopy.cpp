#include "homotopy.h"

extern Random RND;

Homotopy::Homotopy()
{
  fitness = 0.0;
}

Homotopy::Homotopy(const Homotopy& source)
{
  sequence = source.sequence;
  fitness = source.fitness;
}

Homotopy::Homotopy(int n)
{
  for(int i=0; i<n; ++i) {
    Group g;
    sequence.push_back(g);
  }
  compute_fitness();
}

Homotopy::~Homotopy()
{

}

void Homotopy::clear()
{
  sequence.clear();
  fitness = 0.0;
}

void Homotopy::compute_fitness()
{
  unsigned int i;
  double temp,sum = 0.0;
  for(i=0; i<sequence.size(); ++i) {
    temp = std::exp(-std::pow(double(sequence[i].ngenerator)-3.0,2))+ std::pow(std::sin(double(sequence[i].relations.size())),2);
    sum += temp;
  }
  fitness = sum;
}

void Homotopy::mutate()
{
  unsigned int n = RND.irandom(sequence.size());
  Group g;
  sequence[n] = g;
  compute_fitness();
}

std::string Homotopy::write() const 
{
  int i,n = (signed) sequence.size() - 1;
  std::string output = "[";
  for(i=0; i<n; ++i) {
    output += sequence[i].compact_form() + "],[";
  }
  output += sequence[n].compact_form() + "]";
  return output;
}

void Homotopy::serialize(std::ofstream& s) const
{
  int i,n = (signed) sequence.size();
  s.write((char*)(&n),sizeof(int));
  for(i=0; i<n; ++i) {
    sequence[i].serialize(s);
  }
  s.write((char*)(&fitness),sizeof(double));
}

void Homotopy::deserialize(std::ifstream& s)
{
  int i,n;
  Group g;

  clear();

  s.read((char*)(&n),sizeof(int));
  for(i=0; i<n; ++i) {
    g.deserialize(s);
    sequence.push_back(g);
  }
  s.read((char*)(&fitness),sizeof(double));
}

Homotopy& Homotopy::operator =(const Homotopy& source)
{
  if (this == &source) return *this;

  sequence = source.sequence;
  fitness = source.fitness;

  return *this;
}

Homotopy operator +(const Homotopy& h1,const Homotopy& h2)
{
  if (h1.sequence.size() != h2.sequence.size()) {
    std::cerr << "These two homotopy sequences cannot be added: they don't have the same length!" << std::endl;
    std::exit(1);
  }
  unsigned int i,bisection = 1 + RND.irandom(h1.sequence.size()-1);
  Homotopy output;

  for(i=0; i<bisection; ++i) {
    output.sequence.push_back(h1.sequence[i]);
  }
  for(i=bisection; i<h1.sequence.size(); ++i) {
    output.sequence.push_back(h2.sequence[i]);
  }

  return output;
}

void Homotopy::compute(const Nexus* NX)
{
  int i,j,q,r,ngen,n1,n2,n3,e1,e2,nf,vx[3],ntree;
  bool found;
  Word w(0);
  std::vector<int> tree_edges,generator,s1,s2;
  std::vector<Word> relations;
  const int ne = (signed) NX->elements[1].size();
  const int nr = (NX->dimension > 1) ? (signed) NX->elements[2].size() : 0;

  sequence.clear();

  // We need to first compute the zeroth homotopy "group"...
  if (NX->connected()) {
    // Just the trivial group
    sequence.push_back(Group(0));
  }
  else {
    // There isn't really any sort of group structure here, so we'll 
    // just add the "null" group...
    sequence.push_back(Group());
  }

  // First we need to calculate a spanning tree for the 1-skeleton of this complex...
  ntree = NX->spanning_tree(tree_edges);
  for(i=0; i<ne; ++i) {
    s1.push_back(i);
  }
  for(i=0; i<nr; ++i) {
    s2.push_back(i);
  }

  // In princple, this is a group with ngenerators = nedges - ntree/2 and nrelations = ntriangles but we
  // can ignore those 2-simplices all of whose edges are in the spanning tree
  for(i=0; i<ne; ++i) {
    NX->elements[1][s1[i]].get_vertices(vx);

    found = false;
    for(j=0; j<ntree; j+=2) {
      if (vx[0] == tree_edges[j] && vx[1] == tree_edges[j+1]) {
        found = true;
        break;
      }
    }

    if (!(found)) {
      // It's a generator...
      generator.push_back(vx[0]);
      generator.push_back(vx[1]);
    }
  }
  ngen = (signed) generator.size()/2;
  for(i=0; i<nr; ++i) {
    w.clear();
    NX->elements[2][s2[i]].get_vertices(vx);

    nf = 0;

    e1 = vx[0];
    e2 = vx[1];
    found = false;
    n1 = -1;
    for(j=0; j<ngen; ++j) {
      if (e1 == generator[2*j] && e2 == generator[2*j+1]) {
        found = true;
        n1 = j;
        break;
      }
    }
    if (found) nf++;

    e1 = vx[0];
    e2 = vx[2];
    found = false;
    n2 = -1;
    for(j=0; j<ngen; ++j) {
      if (e1 == generator[2*j] && e2 == generator[2*j+1]) {
        found = true;
        n2 = j;
        break;
      }
    }
    if (found) nf++;

    e1 = vx[1];
    e2 = vx[2];
    found = false;
    n3 = -1;
    for(j=0; j<ngen; ++j) {
      if (e1 == generator[2*j] && e2 == generator[2*j+1]) {
        found = true;
        n3 = j;
        break;
      }
    }
    if (found) nf++;

    if (nf == 0) continue;

    if (nf == 1) {
      if (n1 >= 0) w.content.push_back(std::pair<int,int>(n1,1));
      if (n2 >= 0) w.content.push_back(std::pair<int,int>(n2,1));
      if (n3 >= 0) w.content.push_back(std::pair<int,int>(n3,1));
      relations.push_back(w);
      continue;
    }
    // At least two of the edges of this triangle are generators, we need to
    // order them correctly...
    if (n1 >= 0) {
      w.content.push_back(std::pair<int,int>(n1,1));
      q = generator[2*n1+1];

      if (n2 >= 0) {
        r = generator[2*n2];
        if (r == q) {
          w.content.push_back(std::pair<int,int>(n2,1));
          if (n3 >= 0) w.content.push_back(std::pair<int,int>(n3,1));
        }
        else {
          if (n3 >= 0) w.content.push_back(std::pair<int,int>(n3,1));
          w.content.push_back(std::pair<int,int>(n2,-1));
        }
        relations.push_back(w);
        continue;
      }
      r = generator[2*n3];
      if (r == q) {
        w.content.push_back(std::pair<int,int>(n3,1));
      }
      else {
        w.content.push_back(std::pair<int,int>(n3,-1));
      }
    }
    else {
      // So the two non-trivial generators are n2 and n3
      w.content.push_back(std::pair<int,int>(n2,1));
      q = generator[2*n2+1];
      r = generator[2*n3];
      if (q == r) {
        w.content.push_back(std::pair<int,int>(n3,1));
      }
      else {
        w.content.push_back(std::pair<int,int>(n3,-1));
      }
    }
    relations.push_back(w);
  }
  sequence.push_back(Group(ngen,relations));
}

