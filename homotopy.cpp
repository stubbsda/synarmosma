#include "homotopy.h"

using namespace SYNARMOSMA;

extern Random RND;

Homotopy::Homotopy()
{

}

Homotopy::Homotopy(unsigned int n)
{
  if (n == 0) throw std::invalid_argument("The upper limit of the homotopy sequence must be greater than zero!");
  ulimit = n;
  // We set $\pi_0(X)$ to be the trivial group, thereby assuming $X$ to be connected...
  sequence.push_back(Group(1));
  if (ulimit == 1) return;
  // The fundamental group is a random group on six generators with two relations...
  sequence.push_back(Group(6,2));
  if (ulimit == 2) return;
  unsigned int i;
  Group g(10);
  for(i=0; i<ulimit-2; ++i) {
    // The higher-order groups are assumed to be random Abelian groups on ten generators with four relations...
    g.initialize(4);
    g.abelianize();
    sequence.push_back(g);
  }
  compute_fitness();
}

Homotopy::Homotopy(const Homotopy& source)
{
  sequence = source.sequence;
  fitness = source.fitness;
  ulimit = source.ulimit;
}

Homotopy& Homotopy::operator =(const Homotopy& source)
{
  if (this == &source) return *this;

  sequence = source.sequence;
  fitness = source.fitness;
  ulimit = source.ulimit;

  return *this;
}

Homotopy::~Homotopy()
{

}

void Homotopy::clear()
{
  fitness = 0.0;
  ulimit = 0;
  sequence.clear();
}

void Homotopy::compute_fitness()
{
  unsigned int i,ng,nr;
  double sum = 0.0;

  for(i=0; i<ulimit; ++i) {
    ng = sequence[i].get_number_generators();
    nr = sequence[i].get_number_relations();
    sum += std::exp(-std::pow(double(ng - 3),2)) + std::pow(std::sin(double(nr)),2);
  }
  fitness = sum;
}

void Homotopy::mutate(double severity)
{
  if (severity < -std::numeric_limits<double>::epsilon()) throw std::invalid_argument("The mutation severity must be greater than zero!");
  if ((severity - 1.0) > std::numeric_limits<double>::epsilon()) throw std::invalid_argument("The mutation severity must not be greater than one!");
  // If the severity is zero, then do nothing...
  if (severity < std::numeric_limits<double>::epsilon()) return;
  unsigned int n = RND.irandom(1,ulimit);
  unsigned int ng = sequence[n].get_number_generators();
  unsigned int mg = 1 + int(RND.drandom(severity*severity,2.5*severity)*double(ng));
  Group g(mg,mg/2);

  if (n > 1) g = g.abelianize(); 
  sequence[n] = g;
  compute_fitness();
}

std::string Homotopy::write() const
{
  std::string output("[");

  if (ulimit == 0) {
    output += "]";
    return output;
  }
  unsigned int i;

  for(i=0; i<ulimit-1; ++i) {
    output += sequence[i].compact_form() + "],[";
  }
  output += sequence[ulimit-1].compact_form() + "]";
  return output;
}

int Homotopy::serialize(std::ofstream& s) const
{
  unsigned int i;
  int count = 0;

  s.write((char*)(&ulimit),sizeof(int)); count += sizeof(int);
  for(i=0; i<ulimit; ++i) {
    count += sequence[i].serialize(s);
  }

  return count;
}

int Homotopy::deserialize(std::ifstream& s)
{
  unsigned int i;
  int count = 0;
  Group g;

  clear();

  s.read((char*)(&ulimit),sizeof(int)); count += sizeof(int);
  for(i=0; i<ulimit; ++i) {
    count += g.deserialize(s);
    sequence.push_back(g);
  }
  compute_fitness();

  return count;
}

void Homotopy::compute(const Nexus* NX)
{
  int q,r,e1,e2,n1,n2,n3,vx[3];
  unsigned int i,j,nf,ngen,ntree;
  bool found;
  std::vector<int> generator,tree_edges;
  std::vector<unsigned int> s1,s2;
  std::vector<Word> relations;
  const unsigned int ne = NX->get_length(1);
  const unsigned int nr = (NX->get_dimension() > 1) ? NX->get_length(2) : 0;

  clear();
  // Sanity check...
  if (!NX->connected()) return;

  // Since this simplicial complex is connected, its 0-th order homotopy group is 
  // simply the trivial group...
  sequence.push_back(Group(1));
  ulimit += 1;

  // First we need to calculate a spanning tree for the 1-skeleton of this complex...
  ntree = NX->spanning_tree(tree_edges);
  for(i=0; i<ne; ++i) {
    s1.push_back(i);
  }
  for(i=0; i<nr; ++i) {
    s2.push_back(i);
  }

  // In principle, this is a group with ngenerators = nedges - ntree/2 and nrelations = ntriangles but we
  // can ignore those 2-simplices all of whose edges are in the spanning tree
  for(i=0; i<ne; ++i) {
    NX->get_elements(1,s1[i],vx); 

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
  ngen = generator.size()/2;
  Word w;
  for(i=0; i<nr; ++i) {
    w.clear();
    NX->get_elements(2,s2[i],vx); 

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
      if (n1 >= 0) w.append(std::pair<unsigned int,int>(n1,1));
      if (n2 >= 0) w.append(std::pair<unsigned int,int>(n2,1));
      if (n3 >= 0) w.append(std::pair<unsigned int,int>(n3,1));
      relations.push_back(w);
      continue;
    }
    // At least two of the edges of this triangle are generators, we need to
    // order them correctly...
    if (n1 >= 0) {
      w.append(std::pair<unsigned int,int>(n1,1));
      q = generator[2*n1+1];

      if (n2 >= 0) {
        r = generator[2*n2];
        if (r == q) {
          w.append(std::pair<unsigned int,int>(n2,1));
          if (n3 >= 0) w.append(std::pair<unsigned int,int>(n3,1));
        }
        else {
          if (n3 >= 0) w.append(std::pair<unsigned int,int>(n3,1));
          w.append(std::pair<unsigned int,int>(n2,-1));
        }
        relations.push_back(w);
        continue;
      }
      r = generator[2*n3];
      if (r == q) {
        w.append(std::pair<unsigned int,int>(n3,1));
      }
      else {
        w.append(std::pair<unsigned int,int>(n3,-1));
      }
    }
    else {
      // So the two non-trivial generators are n2 and n3
      w.append(std::pair<unsigned int,int>(n2,1));
      q = generator[2*n2+1];
      r = generator[2*n3];
      if (q == r) {
        w.append(std::pair<unsigned int,int>(n3,1));
      }
      else {
        w.append(std::pair<unsigned int,int>(n3,-1));
      }
    }
    relations.push_back(w);
  }
  sequence.push_back(Group(ngen,relations));
  ulimit += 1;
  compute_fitness();
}

namespace SYNARMOSMA {
  std::ostream& operator <<(std::ostream& s,const Homotopy& h)
  {
    s << "[";
    if (h.ulimit == 0) {
      s << "]";
      return s;
    }
    unsigned int i;

    for(i=0; i<h.ulimit-1; ++i) {
      s << h.sequence[i].compact_form() << "],[";
    }
    s << h.sequence[h.ulimit-1].compact_form() << "]";
    return s;
  }

  Homotopy operator +(const Homotopy& h1,const Homotopy& h2)
  {
    if (h1.ulimit != h2.ulimit) throw std::invalid_argument("Homotopy sequences of unequal length cannot be added!");
    unsigned int i,bisection = 1 + RND.irandom(h1.ulimit - 1);
    Homotopy output;

    output.ulimit = h1.ulimit;
    for(i=0; i<bisection; ++i) {
      output.sequence.push_back(h1.sequence[i]);
    }
    for(i=bisection; i<h1.ulimit; ++i) {
      output.sequence.push_back(h2.sequence[i]);
    }
    output.compute_fitness();

    return output;
  }
}
