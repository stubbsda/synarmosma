#include "proposition.h"

using namespace SYNARMOSMA;

extern Random RND;

Proposition::Proposition()
{

}

Proposition::Proposition(const std::set<int>& atoms)
{
  unsigned int nc = 1 + RND.irandom(atoms.size()/Proposition::NP);
  initialize(nc,atoms);
}

Proposition::Proposition(unsigned int nc,const std::set<int>& atoms)
{
  initialize(nc,atoms);
}

Proposition::Proposition(const Proposition& source)
{
  clause = source.clause;
}

Proposition& Proposition::operator =(const Proposition& source)
{
  if (this == &source) return *this;

  clause = source.clause;

  return *this;
}

Proposition::~Proposition()
{

}

void Proposition::initialize(unsigned int nc,const std::set<int>& atoms)
{
  clear();
  if (nc == 0 || atoms.empty()) return;

  int a;
  unsigned int i,j,k;
  double alpha;
  std::set<int> used;
  const unsigned int na = atoms.size();

  for(i=0; i<nc; ++i) {
    a = RND.irandom(atoms);
    clause.push_back(a);
    used.clear();
    used.insert(a);
    if (RND.irandom(2) == 0) {
      clause.push_back(0);
    }
    else {
      clause.push_back(1);
    }
    for(j=1; j<Proposition::NP; ++j) {
      alpha = double(j)/double(2*Proposition::NP);
      if (RND.drandom() < alpha || used.size() == na) {
        for(k=j; k<Proposition::NP; ++k) {
          clause.push_back(-1); clause.push_back(0);
        }
        break;
      }
      a = RND.irandom(atoms,used);
      used.insert(a);
      clause.push_back(a);
      if (RND.irandom(2) == 0) {
        clause.push_back(0);
      }
      else {
        clause.push_back(1);
      }
    }
  }
}

bool Proposition::satisfiable() const
{
  int a;
  unsigned int i,j,k,l,ulimit = 2*Proposition::NP;
  std::unordered_map<int,bool> atom_values;
  std::set<int> atoms,candidates;
  std::set<unsigned int> false_clauses;
  std::set<int>::const_iterator it;
  std::unordered_map<int,bool>::const_iterator qt;

  unsigned int natoms = get_atoms(atoms);
  for(it=atoms.begin(); it!=atoms.end(); ++it) {
    atom_values[*it] = (RND.irandom(2) == 0) ? true : false;
  }
  for(i=0; i<ulimit; ++i) {
    for(j=0; j<3*natoms; ++j) {
      // If the proposition is true with this mapping of the atomic propositions, we're done
      if (evaluate(atom_values,false_clauses)) return true;
      // If it's not true, select a randomly chosen unsatisfied clause of the Boolean formula
      l = RND.irandom(false_clauses);
      // Obtain all the atomic propositions from that clause in the set of candidates
      candidates.clear();
      for(k=0; k<2*Proposition::NP; k+=2) {
        a = clause[2*Proposition::NP*l+k];
        if (a == -1) break;
        candidates.insert(a);
      }
      // Pick one of these atomic propositions randomly from among the candidates and alter its parity
      a = RND.irandom(candidates);
      qt = atom_values.find(a);
      atom_values[a] = !(qt->second);
    }
  }
  return false;
}

void Proposition::mutate()
{
  std::set<int> atoms;
  unsigned int n,nc = get_size();

  n = get_atoms(atoms);
  clear();
  if (n == 0) return;
  initialize(nc,atoms);
}

bool Proposition::evaluate(const std::unordered_map<int,bool>& atoms) const
{
  int a;
  unsigned int i,j;
  bool clause_value,atom_value,output = true; 
  std::unordered_map<int,bool>::const_iterator qt;
  const unsigned int nc = get_size();
 
  for(i=0; i<nc; ++i) {
    // If a single clause is false, the whole expression is false...
    clause_value = false;
    for(j=0; j<2*Proposition::NP; j+=2) {
      a = clause[2*Proposition::NP*i+j];
      if (a == -1) break;
      qt = atoms.find(a);
#ifdef DEBUG
      if (qt == atoms.end()) throw std::invalid_argument("Missing atomic proposition value in Proposition::evaluate method!");
#endif
      atom_value = qt->second;
      if (clause[2*Proposition::NP*i+j+1] == 1) atom_value = !atom_value;
      clause_value = clause_value || atom_value;
    }
    if (!clause_value) {
      output = false;
      break;
    }
  }
  return output;
}

bool Proposition::evaluate(const std::unordered_map<int,bool>& atoms,std::set<unsigned int>& false_clauses) const
{
  int a;
  unsigned int i,j;
  bool clause_value,atom_value,output = true;
  std::unordered_map<int,bool>::const_iterator qt;
  const unsigned int nc = get_size();

  false_clauses.clear();
  for(i=0; i<nc; ++i) {
    // If a single clause is false, the whole expression is false...
    clause_value = false;
    for(j=0; j<2*Proposition::NP; j+=2) {
      a = clause[2*Proposition::NP*i+j];
      if (a == -1) break;
      qt = atoms.find(a);
#ifdef DEBUG
      if (qt == atoms.end()) throw std::invalid_argument("Missing atomic proposition value in Proposition::evaluate method!");
#endif
      atom_value = qt->second;
      if (clause[2*Proposition::NP*i+j+1] == 1) atom_value = !atom_value;
      clause_value = clause_value || atom_value;
    }
    if (!clause_value) {
      false_clauses.insert(i);
      output = false;
    }
  }
  return output;
}

unsigned int Proposition::get_atoms(std::set<int>& atoms) const
{
  int a;
  unsigned int i,j,na = 0,nc = get_size();

  atoms.clear();
  for(i=0; i<nc; ++i) {
    for(j=0; j<2*Proposition::NP; j+=2) {
      a = clause[2*Proposition::NP*i+j];
      if (a == -1) break;
      atoms.insert(a); na++;
    }
  }
  return na;
}

void Proposition::set_atoms(const std::set<int>& atoms)
{
  clear();
  if (atoms.empty()) return;
  if (atoms.size() < Proposition::NP) throw std::invalid_argument("The number of atoms must not be less than the number of atomic propositions!");
  unsigned int nc = 1 + RND.irandom(atoms.size()/Proposition::NP);
  initialize(nc,atoms);
}

int Proposition::serialize(std::ofstream& s) const
{
  int i,count = 0,n = (signed) clause.size();

  s.write((char*)(&n),sizeof(int)); count += sizeof(int);
  for(i=0; i<n; ++i) {
    s.write((char*)(&clause[i]),sizeof(int)); count += sizeof(int);
  }
  return count;
}

int Proposition::deserialize(std::ifstream& s)
{
  int i,j,n,count = 0;

  clear();

  s.read((char*)(&n),sizeof(int)); count += sizeof(int);
  for(i=0; i<n; ++i) {
    s.read((char*)(&j),sizeof(int)); count += sizeof(int);
    clause.push_back(j);
  }
  return count;
}

namespace SYNARMOSMA {
  std::ostream& operator <<(std::ostream& s,const Proposition& p)
  {
    int pvalue,parity;
    unsigned int i,j;
    const unsigned int nc = p.get_size();

    s << "(";
    for(i=0; i<nc; ++i) {
      for(j=0; j<2*(Proposition::NP-1); j+=2) {
        pvalue = p.clause[2*Proposition::NP*i+j];
        parity = p.clause[2*Proposition::NP*i+j+1];
        if (pvalue < 0) continue;
        if (p.clause[2*Proposition::NP*i+j+2] < 0) {
          if (parity == 1) s << "¬";
          s << "p[" << 1 + pvalue << "]";
          continue;
        }
        if (parity == 1) s << "¬";
        s << "p[" << 1 + pvalue << "] ∨ ";
      }
      pvalue = p.clause[2*Proposition::NP*i+2*Proposition::NP-2];
      parity = p.clause[2*Proposition::NP*i+2*Proposition::NP-1];
      if (pvalue >= 0) {
        if (parity == 1) s << "¬";
        s << "p[" << 1 + pvalue << "]";
      }
      if (i < (nc - 1)) {
        s << ") ∧" << std::endl; 
        s << "(";
      }
      else {
        s << ")";
      }
    }
    return s;
  }

  Proposition operator &(const Proposition& P1,const Proposition& P2)
  {
    unsigned int i;
    Proposition output = P1;
    for(i=0; i<P2.clause.size(); ++i) {
      output.clause.push_back(P2.clause[i]);
    }
    return output;
  }
}

