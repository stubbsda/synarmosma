#include "propositional_system.h"

using namespace SYNARMOSMA;

Propositional_System::Propositional_System()
{

}

Propositional_System::Propositional_System(unsigned int na,unsigned int nt)
{
  // nt = number of theorems und na = number of atoms
  if (na == 0 || nt == 0) throw std::invalid_argument("The number of atomic propositions and the number of theorems must be greater than zero!");
  natom = na;
  initialize(nt);
}

Propositional_System::Propositional_System(const Propositional_System& source)
{
  theorems = source.theorems;
  natom = source.natom;
  nuniverse = source.nuniverse;
  truth = source.truth;
}

Propositional_System& Propositional_System::operator =(const Propositional_System& source)
{
  if (this == &source) return *this;

  theorems = source.theorems;
  natom = source.natom;
  nuniverse = source.nuniverse;
  truth = source.truth;

  return *this;
}

Propositional_System::~Propositional_System()
{

}

void Propositional_System::clear()
{
  natom = 0;
  nuniverse = 0;
  truth.clear();
  theorems.clear();
}

void Propositional_System::initialize(unsigned int n)
{
  unsigned int i;
  std::set<int> atoms;

  for(i=0; i<natom; ++i) {
    atoms.insert(i);
  }
  for(i=0; i<n; ++i) {
    theorems.push_back(Proposition(atoms));
  }
  compute_internal_logic();
}

int Propositional_System::serialize(std::ofstream& s) const
{
  unsigned int i,n = theorems.size();
  int count = 0;

  s.write((char*)(&natom),sizeof(int)); count += sizeof(int);
  s.write((char*)(&nuniverse),sizeof(int)); count += sizeof(int);
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
  Proposition p;

  clear();

  s.read((char*)(&natom),sizeof(int)); count += sizeof(int);
  s.read((char*)(&nuniverse),sizeof(int)); count += sizeof(int);
  s.read((char*)(&n),sizeof(int)); count += sizeof(int);
  for(i=0; i<n; ++i) {
    count += p.deserialize(s);
    theorems.push_back(p);
  }
  for(i=0; i<n; ++i) {
    truth.push_back(boost::dynamic_bitset<>(nuniverse));
  }
  compute_internal_logic();

  return count;
}

unsigned int Propositional_System::consistency(unsigned int i,unsigned int j,std::string& type) const
{
  if (i >= theorems.size() || j >= theorems.size()) throw std::invalid_argument("Illegal theorem index in Propositional_System class!");
  boost::dynamic_bitset<> temp(nuniverse);

  boost::to_upper(type);
  if (type == "AND") {
    temp = truth[i] & truth[j];
  }
  else if (type == "OR") {
    temp = truth[i] | truth[j];
  }
  else if (type == "XOR") {
    temp = truth[i] ^ truth[j];
  }
  else {
    throw std::invalid_argument("Unknown Boolean operator!");
  }
  return temp.count();
}

bool Propositional_System::implication(unsigned int n,const std::vector<unsigned int>& axioms) const
{
  if (n >= theorems.size()) throw std::invalid_argument("Illegal theorem index in Propositional_System::implication!");
  unsigned int i;
  boost::dynamic_bitset<> temp(nuniverse);

  temp = truth[n];
  for(i=0; i<axioms.size(); ++i) {
    if (axioms[i] >= theorems.size()) throw std::invalid_argument("Illegal axiom index in Propositional_System::implication!");
    temp |= ~truth[axioms[i]];
  }
  if (temp.count() == nuniverse) return true;
  return false;
}

void Propositional_System::compute_internal_logic()
{
  // The same basic idea as above, fill up the array of edge weights
  // with values but based on propositional logic rather than number
  // theory
  unsigned int i,j,k,in1,p,idiv;
  bool out,row_vector[natom];
  std::unordered_map<int,bool> atom_values;
  Binary_Matrix logical_universe(nuniverse,natom);
  const unsigned int n = theorems.size();

  truth.clear();
  nuniverse = ipow(2,natom);
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
      for(k=0; k<natom; ++k) {
        atom_values[k] = row_vector[k];
      }
      out = theorems[i].evaluate(atom_values);
      if (out) truth[i].set(j);
    }
  }
}

void Propositional_System::compute_implication_graph(Directed_Graph* G) const 
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

unsigned int Propositional_System::add_theorem(const std::set<int>& atoms)
{
  unsigned int n = theorems.size();
  theorems.push_back(Proposition(atoms));
  // Check to see if this proposition uses any new atomic propositions...
  int max_atom = *atoms.rbegin();
  if (max_atom >= natom) {
    natom = max_atom + 1;
    compute_internal_logic();
  }
  else {
    truth.push_back(boost::dynamic_bitset<>(nuniverse));
  }
  return n;
}

bool Propositional_System::drop_theorem(unsigned int n)
{
  if (n >= theorems.size()) return false;
  theorems.erase(theorems.begin() + n);
  truth.erase(truth.begin() + n);
  return true;
}