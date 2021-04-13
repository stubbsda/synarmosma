#include "poset.h"

using namespace SYNARMOSMA;

Poset::Poset()
{

}

Poset::Poset(int n)
{
  // We begin by assuming that every pair of points is
  // disparate, so the hash map is empty
  if (n < 1) throw std::invalid_argument("A poset must have at least one element!");
  N = n;
}

Poset::Poset(int n,double lambda)
{
  if (n < 1) throw std::invalid_argument("A poset must have at least one element!");
  N = n;
  construct_ordering(lambda);
}

Poset::~Poset()
{

}

Poset::Poset(const Poset& source)
{
  N = source.N;
  order = source.order;
}

Poset& Poset::operator =(const Poset& source) 
{
  if (this == &source) return *this;

  N = source.N;
  order = source.order;

  return *this;
}

void Poset::clear()
{
  N = 0;
  order.clear();
}

int Poset::serialize(std::ofstream& s) const
{
  int i,j;
  int count = 0;
  Relation rho;

  s.write((char*)(&N),sizeof(int)); count += sizeof(int);
  for(i=0; i<N; ++i) {
    for(j=1+i; j<N; ++j) {
      rho = get_order(i,j);
      s.write((char*)(&rho),sizeof(int)); count += sizeof(int);
    }
  }

  return count;
}

int Poset::deserialize(std::ifstream& s)
{
  int i,j;
  int count = 0;
  Relation rho;

  clear();

  s.read((char*)(&N),sizeof(int)); count += sizeof(int);
  for(i=0; i<N; ++i) {
    for(j=1+i; j<N; ++j) {
      s.read((char*)(&rho),sizeof(int)); count += sizeof(int);
      if (rho == Relation::before) {
        set_order(i,j);
      }
      else if (rho == Relation::after) {
        set_order(j,i);
      }
    }
  }
  if (!consistent()) throw std::runtime_error("The poset created using the deserialize method is inconsistent!");

  return count;
}

bool Poset::invert_order(int u,int v)
{
  if (u == v) return false;
  Relation rho = get_order(u,v);
  if (rho == Relation::disparate) return false;
  if (rho == Relation::before) {
    // Inverting order of u < v
    set_order(v,u);
  }
  else {
    // Inverting order of v < u
    set_order(u,v);
  }
  return true;
}

bool Poset::unset_order(int u,int v)
{
  if (u == v) return false;

  Relation rho = get_order(u,v);
  if (rho == Relation::disparate) return false;
  if (rho == Relation::before) {
    // Removing order of u < v
    order.erase(std::pair<int,int>(u,v));
    for(int i=0; i<N; ++i) {
      if (i == u || i == v) continue;
      if (get_order(u,i) == Relation::before && get_order(i,v) == Relation::before) unset_order(i,v);
    }
  }
  else {
    // Removing order of v < u
    order.erase(std::pair<int,int>(v,u));
    for(int i=0; i<N; ++i) {
      if (i == u || i == v) continue;
      if (get_order(v,i) == Relation::before && get_order(i,u) == Relation::before) unset_order(i,u);
    }
  }
  return true;
}

bool Poset::set_order(int u,int v)
{
  if (u == v) return false;
  // Check to see if we already have u < v or v < u in this poset
  Relation rho = get_order(u,v);
  if (rho == Relation::before) return false;

  if (rho == Relation::after) order.erase(std::pair<int,int>(v,u));
  order[std::pair<int,int>(u,v)] = true;
  // Keep the ordering transitive!
  for(int i=0; i<N; ++i) {
    if (i == u || i == v) continue;
    if (get_order(v,i) == Relation::before) set_order(u,i);
    if (get_order(i,u) == Relation::before) set_order(i,v);
  }
  return true;
}

Relation Poset::get_order(int u,int v) const
{
  if (u == v) return Relation::before;
  boost::unordered_map<std::pair<int,int>,bool>::const_iterator qt;
  qt = order.find(std::pair<int,int>(u,v));
  if (qt != order.end()) return Relation::before;
  qt = order.find(std::pair<int,int>(v,u));
  if (qt != order.end()) return Relation::after;
  return Relation::disparate;
}

int Poset::build_chain(std::vector<int>& chain,int length) const
{
  if (chain.empty()) throw std::invalid_argument("The chain vector in Poset::build_chain must not be empty!");

  int l = (signed) chain.size(),output = 0;
  if (l == length) {
    output = 1;
  }
  else {
    int i;
    std::vector<int> nchain = chain;

    for(i=0; i<N; ++i) {
      if (i == chain[l]) continue;
      if (get_order(chain[l],i) == Relation::before) {
        nchain.push_back(i);
        output += build_chain(nchain,length);
        nchain = chain;
      }
    }
  }
  return output;
}

int Poset::chain_number(int length) const
{
  // Compute the number of chains of a given length in this poset 
  if (length < 0) throw std::invalid_argument("The length of the chain in the Poset class must be non-negative!");
  if (length == 0 || length == 1) return 0;
  int i,j,nchain = 0;
  std::set<int> sigma;

  if (length == 2) {
    // The easiest case
    for(i=0; i<N; ++i) {
      for(j=0; j<N; ++j) {
        if (i == j) continue;
        if (get_order(i,j) == Relation::before) nchain++;
      }
    }
  }
  else if (length == 3) {
    // Fairly easy as well
    for(i=0; i<N; ++i) {
      for(j=0; j<N; ++j) {
        if (i == j) continue;
        if (get_order(i,j) != Relation::before) continue;
        compute_width(i,j,sigma);
        nchain += (signed) sigma.size();
      }
    }
  }
  else {
    // The general case, rather complicated so we use recursion to do it
    std::vector<int> chain;
   
    for(i=0; i<N; ++i) {
      chain.push_back(i);
      nchain += build_chain(chain,length);
      chain.clear();
    }   
  }
  return nchain;
}

void Poset::compute_width(int u,int v,std::set<int>& slice) const
{
  if (u == v) throw std::invalid_argument("The two elements in Poset::compute_width must be distinct!");

  slice.clear();

  for(int i=0; i<N; ++i) {
    if (i == u || i == v) continue;
    if (get_order(u,i) != Relation::before) continue;
    if (get_order(i,v) != Relation::before) continue;
    slice.insert(i);
  }
}

void Poset::power_set(int n)
{
  // A method that builds a poset based on the inclusion relation for 
  // the power set on n elements.
  int k = 0;
  bool inclusion;
  int i,j,stack[1+n];
  std::set<int> S;
  std::set<int>::const_iterator it;
  std::vector<std::set<int> > powerset;

  clear();
  N = ipow(2,n);

  stack[0] = 0; 
  powerset.push_back(S);
  while(true) {
    if (stack[k] < n) {
      stack[k+1] = stack[k] + 1;
      k++;
    }
    else {
      stack[k-1] = stack[k-1] + 1;
      k--;
    }

    if (k == 0) break;
    S.clear();
    for(i=1; i<=k; ++i) {
      S.insert(stack[i]);
    }
    powerset.push_back(S);
  }

  for(i=0; i<N; ++i) {
    for(j=1+i; j<N; ++j) {
      inclusion = true;
      for(it=powerset[i].begin(); it!=powerset[i].end(); ++it) {
        if (powerset[j].count(*it) == 0) {
          inclusion = false;
          break;
        }
      }
      if (inclusion) {
        set_order(i,j);
        continue;
      }
      inclusion = true;
      for(it=powerset[j].begin(); it!=powerset[j].end(); ++it) {
        if (powerset[i].count(*it) == 0) {
          inclusion = false;
          break;
        }
      }
      if (inclusion) set_order(j,i);
    }
  }   
}

double Poset::totality() const
{
  // A method to calculate the percentage of pair-wise relationships in this 
  // poset, i.e. how close is this partial order to being a total order?
  int i,j,norder = 0;
  Relation rho;
  const int ntotal = N*(N-1)/2;

  for(i=0; i<N; ++i) {
    for(j=1+i; j<N; ++j) {
      rho = get_order(i,j);
      if (rho != Relation::disparate) norder++;
    }
  }
  return double(norder)/double(ntotal);
}

bool Poset::covered(int u,int v) const
{
  // A method to determine if u is covered by v
  if (u == v) return false;
  if (get_order(u,v) != Relation::before) return false;
  std::set<int> sigma;
  compute_width(u,v,sigma);
  return sigma.empty(); 
}

bool Poset::compute_covering_graph(Graph* G) const
{
  int i,j;

  G->clear();

  for(i=0; i<N; ++i) {
    G->add_vertex();
  }
  for(i=0; i<N; ++i) {
    for(j=1+i; j<N; ++j) {
      if (covered(i,j) || covered(j,i)) G->add_edge(i,j);
    }
  }
  if (G->connected()) return G->planar();
  return false;
}

bool Poset::sink(int n) const
{
  // This is an element whose posteriority is null
  std::set<int> S;
  compute_posteriority(n,S);
  return S.empty();
}

bool Poset::source(int n) const
{
  // This is an element whose anteriority is null
  std::set<int> S;
  compute_anteriority(n,S);
  return S.empty();
}

void Poset::compute_anteriority(int n,std::set<int>& output) const
{
  int i;
  output.clear();
  for(i=0; i<N; ++i) {
    if (i == n) continue;
    if (get_order(i,n) == Relation::before) output.insert(i);
  }
}

void Poset::compute_posteriority(int n,std::set<int>& output) const
{
  int i;
  output.clear();
  for(i=0; i<N; ++i) {
    if (i == n) continue;
    if (get_order(n,i) == Relation::before) output.insert(i);
  }
}

bool Poset::consistent() const
{
  // We need to make sure the order relation satisfies the axioms of a 
  // partial order, i.e. reflexive, anti-symmetric and transitive. 
  int i,j,k;
  boost::unordered_map<std::pair<int,int>,bool>::const_iterator qt;

  for(i=0; i<N; ++i) {
    // Find every element that is after this one and make sure that all the 
    // elements after them are also after "i"
    for(j=0; j<N; ++j) {
      if (i == j) continue;
      if (get_order(i,j) != Relation::before) continue;
      if (get_order(j,i) != Relation::after) return false;
      for(k=0; k<N; ++k) {
        if (k == j || k == i) continue;
        // So if i < j and j < k, then it must be that i < k
        if (get_order(j,k) == Relation::before) {
          if (get_order(i,k) != Relation::before) return false;
        }
      }
    }
  }
  return true;   
}

void Poset::write_incastrature(const std::string& filename) const
{
  // A method that generates the Hasse diagram corresponding to the
  // complex's simplicial structure, with the diagram stored in the
  // PDF "hasse.pdf"; the method assumes that the Graphviz library
  // has been installed on the system.
  int i,j;

  std::ofstream s(filename,std::ios::trunc);

  s << "digraph G {" << std::endl;
  // First all the elements in the poset...
  for(i=0; i<N; ++i) {
    s << "  \"" << 1+i << "\";" << std::endl;
  }
  // Now the directed edges induced by the poset's ordering...
  for(i=0; i<N; ++i) {
    for(j=0; j<N; ++j) {
      if (i == j) continue;
      if (get_order(i,j) == Relation::before) {
        s << "  \"" << 1+i << "\" -> \"" << 1+j << "\";" << std::endl;
      }
    }
  }
  s << "}" << std::endl;
  s.close();
}

void Poset::construct_ordering(double lambda)
{
  // A method to impose a random order on the poset 
  if (N < 2) throw std::runtime_error("Cannot construct order on a poset with fewer than two elements!");
  if (lambda < 0.0 || lambda > 1.0) throw std::invalid_argument("The desired poset totality must be between 0 and 1!");
  int u,v;
  double p = 0.0;
  Random RND;

  do {
    u = RND.irandom(N);
    v = RND.irandom(N);
    if (set_order(u,v)) p = totality();
  } while(p < lambda);

  if (!consistent()) std::runtime_error("The poset created in Poset::construct_ordering does not have a consistent ordering!");
}

