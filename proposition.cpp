#include "proposition.h"

using namespace SYNARMOSMA;

extern Random RND;

Proposition::Proposition()
{

}

Proposition::Proposition(const std::set<int>& atoms)
{
  int nc = 1 + RND.irandom(atoms.size()/Proposition::NP);
  initialize(nc,atoms);
}

Proposition::Proposition(int nc,const std::set<int>& atoms)
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
  int p;
  unsigned int i,j,k;
  double alpha;
  std::set<int> used;
  const unsigned int na = atoms.size();

  clause.clear();

  for(i=0; i<nc; ++i) {
    p = RND.irandom(atoms);
    clause.push_back(p);
    used.clear();
    used.insert(p);
    if (RND.drandom() < 0.5) {
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
      p = RND.irandom(atoms,used);
      used.insert(p);
      clause.push_back(p);
      if (RND.drandom() < 0.5) {
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
  int i,l,a,p,K = 5;
  unsigned int j,k;
  std::vector<int> bvalues,fclause;
  std::set<int> atoms,ca;
  std::set<int>::const_iterator it;

  get_atoms(atoms);
  for(it=atoms.begin(); it!=atoms.end(); ++it) {
    bvalues.push_back(*it);
    bvalues.push_back(RND.irandom(2));
  }
  const unsigned int natoms = atoms.size();
  for(i=0; i<K; ++i) {
    for(j=0; j<3*natoms; ++j) {
      // Now grab an unsatisfied clause of the Boolean
      // formula
      if (evaluate(bvalues,fclause)) return true;
      l = RND.irandom(fclause);
      for(k=0; k<2*Proposition::NP; k+=2) {
        p = clause[2*Proposition::NP*l+k];
        if (p >= 0) ca.insert(p);
      }
      a = RND.irandom(ca);
      for(k=0; k<2*natoms; k+=2) {
        if (bvalues[k] == a) {
          bvalues[k+1] = (bvalues[k+1] + 1)%2;
          break;
        }
      }
    }
  }
  return false;
}

void Proposition::mutate()
{
  std::set<int> atoms;
  unsigned int nc = get_size();

  get_atoms(atoms);
  clear();
  if (atoms.empty()) return;
  initialize(nc,atoms);
}

void Proposition::mutate(const std::set<int>& atoms)
{
  clear();
  if (atoms.empty()) return;
  unsigned int nc = 1 + RND.irandom(atoms.size()/Proposition::NP); 
  initialize(nc,atoms);
}

bool Proposition::evaluate(const std::vector<std::pair<int,bool> >& atom_values) const
{
  int a;
  unsigned int i,j,k;
  bool clause_value,atom_value,atom_found,output = true; 
  const unsigned int in_size = atom_values.size();
  const unsigned int nc = get_size();
 
  for(i=0; i<nc; ++i) {
    // If a single clause is false, the whole expression is false...
    clause_value = false;
    for(j=0; j<2*Proposition::NP; j+=2) {
      a = clause[2*Proposition::NP*i+j];
      if (a == -1) break;
      atom_found = false;
      for(k=0; k<in_size; ++k) {
        if (atom_values[k].first == a) {
          atom_value = atom_values[k].second;
          atom_found = true;
          break;
        }
      }
      if (!atom_found) throw std::invalid_argument("Missing atomic proposition value in Proposition::evaluate method!");
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

bool Proposition::evaluate(const std::vector<int>& atoms,std::vector<int>& fclause) const
{
  int a;
  unsigned int i,j,k;
  bool cvalue,avalue = false,output = true;
  const unsigned int na = atoms.size();
  const unsigned int nc = get_size();


  fclause.clear();
  for(i=0; i<nc; ++i) {
    // If a single clause is false, the whole expression is false...
    cvalue = false;
    for(j=0; j<2*Proposition::NP; j+=2) {
      a = clause[2*Proposition::NP*i+j];
      if (a == -1) break;
      for(k=0; k<na; k+=2) {
        if (atoms[k] == a) {
          avalue = (atoms[k+1] == 0) ? false : true;
          break;
        }
      }
      if (clause[2*Proposition::NP*i+j+1] == 1) avalue = !avalue;
      cvalue = cvalue || avalue;
    }
    if (!cvalue) {
      fclause.push_back(i);
      output = false;
    }
  }
  return output;
}

void Proposition::get_atoms(std::set<int>& atoms) const
{
  unsigned int i,l = clause.size();

  atoms.clear();
  for(i=0; i<l; i+=2) {
    if (clause[i] < 0) continue;
    atoms.insert(clause[i]);
  }
}

void Proposition::set_atoms(const std::set<int>& atoms)
{
  int nc = 1 + RND.irandom(atoms.size()/Proposition::NP);
  initialize(nc,atoms);
}

int Proposition::serialize(std::ofstream& s) const
{
  int i,count = 0,n = clause.size();

  s.write((char*)(&n),sizeof(int)); count += sizeof(int);
  for(i=0; i<n; ++i) {
    s.write((char*)(&clause[i]),sizeof(int)); count += sizeof(int);
  }
  return count;
}

int Proposition::deserialize(std::ifstream& s)
{
  int i,n,l,count = 0;

  clear();

  s.read((char*)(&n),sizeof(int)); count += sizeof(int);
  for(i=0; i<n; ++i) {
    s.read((char*)(&l),sizeof(int)); count += sizeof(int);
    clause.push_back(l);
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
          if (parity == 1) s << "!";
          s << "p[" << 1 + pvalue << "]";
          continue;
        }
        if (parity == 1) s << "!";
        s << "p[" << 1 + pvalue << "] | ";
      }
      pvalue = p.clause[2*Proposition::NP*i+2*Proposition::NP-2];
      parity = p.clause[2*Proposition::NP*i+2*Proposition::NP-1];
      if (pvalue >= 0) {
        if (parity == 1) s << "!";
        s << "p[" << 1 + pvalue << "]";
      }
      if (i < (nc - 1)) {
        s << ") &" << std::endl; 
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

