#include "propositional_system.h"

using namespace SYNARMOSMA;

Propositional_System::Propositional_System()
{

}

Propositional_System::Propositional_System(const std::set<int>& a,unsigned int nt)
{
  if (a.empty()) throw std::invalid_argument("The number of atomic propositions must be greater than zero!");
  if (nt == 0) throw std::invalid_argument("The number of theorems must be greater than zero!");
  atoms = a;
  initialize(nt);
}

Propositional_System::Propositional_System(const Propositional_System& source)
{
  theorems = source.theorems;
  atoms = source.atoms;
  nuniverse = source.nuniverse;
}

Propositional_System& Propositional_System::operator =(const Propositional_System& source)
{
  if (this == &source) return *this;

  theorems = source.theorems;
  atoms = source.atoms;
  nuniverse = source.nuniverse;

  return *this;
}

Propositional_System::~Propositional_System()
{

}

void Propositional_System::clear()
{
  nuniverse = 0;
  atoms.clear();
  theorems.clear();
}

void Propositional_System::initialize(unsigned int n)
{
  unsigned int i;

  for(i=0; i<n; ++i) {
    theorems.push_back(Proposition(atoms)); 
  }
}

int Propositional_System::serialize(std::ofstream& s) const
{
  unsigned int i,n = theorems.size();
  int count = 0;

  s.write((char*)(&n),sizeof(int)); count += sizeof(int);
  for(i=0; i<n; ++i) {
    count += theorems[i].serialize(s);
  }
  return count;
}

int Propositional_System::deserialize(std::ifstream& s)
{
  unsigned int i,n;
  int count = 0;
  std::set<int> a;
  std::set<int>::const_iterator it;
  Proposition p;

  clear();

  s.read((char*)(&n),sizeof(int)); count += sizeof(int);
  for(i=0; i<n; ++i) {
    count += p.deserialize(s);
    theorems.push_back(p);
    theorems[i].get_atoms(a);
    for(it=a.begin(); it!=a.end(); ++it) {
      atoms.insert(*it);
    }
  }
  nuniverse = UINT64(ipow(2,atoms.size()));

  return count;
}

unsigned int Propositional_System::consistency(unsigned int i,unsigned int j,const std::string& type,const std::vector<boost::dynamic_bitset<> >& truth) const
{
  if (i >= theorems.size() || j >= theorems.size()) throw std::invalid_argument("Illegal theorem index in the Propositional_System::consistency method!");
  boost::dynamic_bitset<> temp(nuniverse);

  std::string utype = boost::to_upper_copy(type);

  if (utype == "AND") {
    temp = truth[i] & truth[j];
  }
  else if (utype == "OR") {
    temp = truth[i] | truth[j];
  }
  else if (utype == "XOR") {
    temp = truth[i] ^ truth[j];
  }
  else {
    throw std::invalid_argument("Unrecognized logical operator in the Propositional_System::consistency method!");
  }
  return temp.count();
}

bool Propositional_System::implication(unsigned int n,const std::vector<unsigned int>& axioms,const std::vector<boost::dynamic_bitset<> >& truth) const
{
  if (n >= theorems.size()) throw std::invalid_argument("Illegal theorem index in the Propositional_System::implication method!");
  unsigned int i;
  boost::dynamic_bitset<> temp(nuniverse);

  temp = truth[n];
  for(i=0; i<axioms.size(); ++i) {
    if (axioms[i] >= theorems.size()) throw std::invalid_argument("Illegal axiom index in the Propositional_System::implication method!");
    temp |= ~truth[axioms[i]];
  }
  if (temp.count() == nuniverse) return true;
  return false;
}

void Propositional_System::compute_internal_logic(std::vector<boost::dynamic_bitset<> >& truth)
{
  // The same basic idea as above, fill up the array of edge weights
  // with values but based on propositional logic rather than number
  // theory
  const unsigned int natom = atoms.size();
  if (natom == 0) return;
  unsigned int i,j,k,in1,p,idiv;
  bool out,row_vector[natom];
  std::set<int>::const_iterator it;
  std::unordered_map<int,bool> atom_values;
  Binary_Matrix logical_universe(nuniverse,natom);
  const unsigned int n = theorems.size();

  truth.clear();
  nuniverse = UINT64(ipow(2,natom));
  for(i=0; i<theorems.size(); ++i) {
    truth.push_back(boost::dynamic_bitset<>(nuniverse));
  }

  // Loop through all vertices, if the topology has been modified, then
  // compute its truth value (do the set of propositions of its neighbours,
  // taken as axioms, imply this vertex's proposition as a logical
  // conclusion), which will be used to compute the edge weights
  for(i=0; i<nuniverse; ++i) {
    in1 = i;
    for(j=0; j<natom; ++j) {
      p = ipow(2,natom-(j-1));
      idiv = in1/p;
      if (!(idiv % 2)) logical_universe.set(i,j);
      in1 -= p*idiv;
    }
  }

  for(i=0; i<n; ++i) {
    truth[i].reset();
    for(j=0; j<nuniverse; ++j) {
      logical_universe.get_row(j,row_vector);
      k = 0;
      for(it=atoms.begin(); it!=atoms.end(); ++it) {
        atom_values[*it] = row_vector[k]; ++k;
      }
      out = theorems[i].evaluate(atom_values);
      if (out) truth[i].set(j);
    }
  }
}

void Propositional_System::compute_implication_graph(Directed_Graph* G,const std::vector<boost::dynamic_bitset<> >& truth) const
{
  // This method must determine which propositions imply other propositions among
  // our collection, "theorems".
  int i,j;
  boost::dynamic_bitset<> temp(nuniverse),p1(nuniverse);

  const int nnode = (signed) theorems.size();

  G->clear();

  for(i=0; i<nnode; ++i) {
    p1 = truth[i];
    for(j=0; j<nnode; ++j) {
      if (i == j) continue;
      // p(i) implies p(j)? That means that (!p(j) v p(i)) is always true
      temp = ~truth[j] | p1;
      if (temp.count() == nuniverse) G->add_edge(i,j,Relation::before);
    }
  }
}

bool Propositional_System::compute_propositional_lattice(Lattice* L) const
{
  int i,j,zero,unity;
  bool subset;
  std::set<int> U,V;
  std::set<int>::const_iterator it;
  const int nnode = (signed) theorems.size();

  L->clear();

  for(i=0; i<nnode; ++i) {
    theorems[i].get_atoms(U);
    if (U.empty()) throw std::runtime_error("Empty theorem in the Propositional_System::compute_propositional_lattice method!");    
  }

  for(i=0; i<nnode; ++i) {
    L->add_element();
  }
  zero = L->add_element();
  unity = L->add_element();

  for(i=0; i<nnode; ++i) {
    theorems[i].get_atoms(U);
    for(j=1+i; j<nnode; ++j) {
      // Check if one of these sets of propositional atoms is a 
      // subset of the other
      theorems[j].get_atoms(V);
      if (U.size() == V.size()) continue;
      if (U.size() < V.size()) {
        subset = true;
        for(it=U.begin(); it!=U.end(); ++it) {
          if (V.count(*it) == 0) {
            subset = false;
            break;
          }
        }
        if (subset) L->set_order(i,j);
      }
      else {
        subset = true;
        for(it=V.begin(); it!=V.end(); ++it) {
          if (U.count(*it) == 0) {
            subset = false;
            break;
          }
        }
        if (subset) L->set_order(j,i);
      }
    }
  }

  for(i=0; i<nnode; ++i) {
    // The zero element is the empty set containing no atomic propositions...
    if (L->source(i)) L->set_order(zero,i);
    // The unity element is the set of all atomic propositions used in the theorems 
    // of this propositional system...
    if (L->sink(i)) L->set_order(i,unity);
  }

  return L->consistent();
}

unsigned int Propositional_System::add_theorem(const std::set<int>& a)
{
  std::set<int>::const_iterator it;
  unsigned int n = theorems.size();

  theorems.push_back(Proposition(a));

  // Check to see if this proposition uses any new atomic propositions...
  for(it=a.begin(); it!=a.end(); ++it) {
    atoms.insert(*it);
  }
  nuniverse = UINT64(ipow(2,atoms.size()));

  return n;
}
