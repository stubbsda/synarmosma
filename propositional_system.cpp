#include "propositional_system.h"

using namespace SYNARMOSMA;

Propositional_System::Propositional_System()
{
  set_default_values();
  initialize(100);
}

Propositional_System::Propositional_System(unsigned int n)
{
  set_default_values();
  initialize(n);
}

Propositional_System::Propositional_System(unsigned int n,const char* filename)
{
  unsigned int i,eq_point;
  std::string line,name,value;

  set_default_values();

  // Open the parameter file
  std::ifstream s(filename,std::ios::in);
  if (s.is_open() == false) {
    // If the file isn't there, print an error message and
    // exit cleanly
    std::cerr << "The file " << filename << " cannot be found!" << std::endl;
    std::exit(1);
  }

  while(std::getline(s,line)) {
    // If it's an empty line, continue
    if (line.empty()) continue;
    // If the line begins with a #, ignore it
    if (line[0] == '#') continue;
    // Find the position of the equals sign
    eq_point = 0;
    for(i=0; i<line.length(); ++i) {
      if (line[i] == '=') {
        eq_point = i;
        break;
      }
    }
    // If there's no equals sign in this line, continue
    if (eq_point < 1) continue;
    name = line.substr(0,eq_point-1);
    trim(name);
    value = line.substr(eq_point+1,line.length());
    trim(value);
    // Now that we have the parameter name, see if it matches
    // any of the known parameters. If so, read in the value and
    // assign it
    if (name == "naxiom") {
      natom = boost::lexical_cast<int>(value);
    }
  }
  initialize(n);
}

Propositional_System::Propositional_System(unsigned int n,unsigned int m)
{
  // n = number of axioms und m = number of atoms
  set_default_values();
  natom = m;
  initialize(n);
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
  theorems.clear();
  natom = 0;
  nuniverse = 0;
  truth.clear();
}

void Propositional_System::set_default_values()
{
  natom = 10;
}

void Propositional_System::initialize(unsigned int n)
{
  unsigned int i;
  std::set<int> atoms;

  for(i=0; i<n; ++i) {
    atoms.insert(i);
  }
  for(i=0; i<n; ++i) {
    theorems.push_back(Proposition(atoms));
  }
  nuniverse = ipow(2,natom);
  for(i=0; i<theorems.size(); ++i) {
    truth.push_back(boost::dynamic_bitset<>(nuniverse));
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

unsigned int Propositional_System::bit_count(unsigned int i) const
{
  return truth[i].count();
}

unsigned int Propositional_System::consistency(unsigned int i,unsigned int j,const std::string& op) const
{
  boost::dynamic_bitset<> temp(nuniverse);
  if (op == "and") {
    temp = truth[i] & truth[j];
  }
  else if (op == "xor") {
    temp = truth[i] ^ truth[j];
  }
  else {
    temp = truth[i] | truth[j];
  }
  return temp.count();
}

bool Propositional_System::implication(unsigned int in1,const std::vector<unsigned int>& axioms) const
{
  unsigned int i;
  boost::dynamic_bitset<> temp(nuniverse);

  temp = truth[in1];
  for(i=0; i<axioms.size(); ++i) {
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
  std::vector<std::pair<int,bool> > rvector;
  Binary_Matrix logical_universe(nuniverse,natom);
  const unsigned int n = theorems.size();

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
      rvector.clear();
      for(k=0; k<natom; ++k) {
        rvector.push_back(std::pair<int,bool>(k,row_vector[k]));
      }
      out = theorems[i].evaluate(rvector);
      if (out) truth[i].set(j);
    }
  }
}

void Propositional_System::compute_implication_graph(Directed_Graph* G) const 
{
  // This method must determine which propositions imply other propositions among
  // our collection, "theorems".
  unsigned int i,j;
  boost::dynamic_bitset<> temp(nuniverse),p1(nuniverse);
  const unsigned int nnode = theorems.size();

  G->clear();

  for(i=0; i<nnode; ++i) {
    p1 = truth[i];
    for(j=0; j<nnode; ++j) {
      if (i == j) continue;
      // p(i) implies p(j)? That means that (!p(j) v p(i)) is always true
      temp = ~truth[j] | p1;
      if (temp.count() == nuniverse) G->add_edge(i,j,OUTGOING);
    }
  }
}



