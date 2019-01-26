#include "lattice.h"

using namespace SYNARMOSMA;

extern Random RND;

Lattice::Lattice() : Poset()
{

}

Lattice::Lattice(int n) : Poset(n)
{
  if (n < 2) throw std::invalid_argument("A lattice must have at least two elements!");
  initialize();
}

Lattice::~Lattice()
{

}

Lattice::Lattice(const Lattice& source)
{
  N = source.N;
  order = source.order;
  atomic = source.atomic;
  atoms = source.atoms;
  null = source.null;
  unity = source.unity;
}

Lattice& Lattice::operator =(const Lattice& source) 
{
  if (this == &source) return *this;

  N = source.N;
  order = source.order;
  atomic = source.atomic;
  atoms = source.atoms;
  null = source.null;
  unity = source.unity;

  return *this;
}

void Lattice::clear()
{
  N = 0;
  order.clear();
  atomic = false;
  null = 0;
  unity = 0;
  atoms.clear();
}

int Lattice::serialize(std::ofstream& s) const
{
  int i,j,count = 0;
  std::set<int>::const_iterator it;
  Relation rho;

  s.write((char*)(&N),sizeof(int)); count += sizeof(int);
  s.write((char*)(&atomic),sizeof(bool)); count += sizeof(bool);
  s.write((char*)(&null),sizeof(int)); count += sizeof(int);
  s.write((char*)(&unity),sizeof(int)); count += sizeof(int);
  for(i=0; i<N; ++i) {
    for(j=1+i; j<N; ++j) {
      rho = get_order(i,j);
      s.write((char*)(&rho),sizeof(int)); count += sizeof(int);
    }
  }
  j = atoms.size();
  s.write((char*)(&j),sizeof(int)); count += sizeof(int);
  for(it=atoms.begin(); it!=atoms.end(); ++it) {
    i = *it;
    s.write((char*)(&i),sizeof(int)); count += sizeof(int);
  }

  return count;
}

int Lattice::deserialize(std::ifstream& s)
{
  int i,j,k,count = 0;
  Relation rho;

  clear();

  s.read((char*)(&N),sizeof(int)); count += sizeof(int);
  s.read((char*)(&atomic),sizeof(bool)); count += sizeof(bool);
  s.read((char*)(&null),sizeof(int)); count += sizeof(int);
  s.read((char*)(&unity),sizeof(int)); count += sizeof(int);
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
  s.read((char*)(&j),sizeof(int)); count += sizeof(int);
  for(i=0; i<j; ++i) {
    s.read((char*)(&k),sizeof(int)); count += sizeof(int);
    atoms.insert(k);
  }

  if (!consistent()) throw std::runtime_error("The lattice created using the deserialize method is inconsistent!");

  return count;
}


void Lattice::initialize()
{
  if (N < 2) throw std::runtime_error("A lattice must have at least two elements!");

  int i,j,n1,n2,ndelta,delta = 2*N*(N-1);

  do {
    n1 = RND.irandom(N);
    n2 = RND.irandom(N);
    if (n1 == n2) continue;
    if (get_order(n1,n2) != Relation::disparate) continue;
    set_order(n1,n2);
    
    ndelta = 0;
    for(i=0; i<N; ++i) {
      for(j=0; j<N; ++j) {
        if (i == j) continue;
        if (meet(i,j) == N) ndelta++;
        if (join(i,j) == N) ndelta++;
      }
    }
    if (ndelta > delta) {
      unset_order(n1,n2);
      continue;
    }
    delta = ndelta;
  } while(delta > 0);

  if (!consistent()) throw std::runtime_error("The lattice created by the initialize method is inconsistent!");

  compute_bounds();
  compute_atoms();
}

bool Lattice::consistent() const
{
  // In a lattice every pair of elements must have a meet and join, so...
  int i,j;

  for(i=0; i<N; ++i) {
    for(j=0; j<N; ++j) {
      if (i == j) continue;
      if (meet(i,j) == N) return false;
      if (join(i,j) == N) return false;
    }
  }
  return true;
}

double Lattice::atomicity() const
{
  if (atomic) return 1.0;

  int i,j,npotential = 0,non_atomic = 0;
  bool found;
  std::set<int>::const_iterator it;

  for(i=0; i<N; ++i) {
    if (i == null) continue;
    if (atoms.count(i) > 0) continue;
    npotential++;
    found = false;
    for(it=atoms.begin(); it!=atoms.end(); ++it) {
      j = *it;
      if (get_order(j,i) == Relation::before) {
        found = true;
        break;
      }
    }
    if (!found) non_atomic++;  
  }
  // Since the lattice is non-atomic, we know that npotential > 0
  return double(non_atomic)/double(npotential);  
}

void Lattice::compute_bounds()
{
  // This method computes the null and identity elements for this lattice using a 
  // brute force technique
  int i,j;
  bool found,nfound;
  
  // We begin by finding the 0 element, i.e. the element of the lattice such that 0 ~ x for all x
  nfound = false;
  for(i=0; i<N; ++i) {
    found = true;
    for(j=0; j<N; ++j) {
      if (i == j) continue;
      if (get_order(i,j) != Relation::before) {
        found = false;
        break;
      }
    }
    if (found) {
      null = i;
      nfound = true;
      break;
    }
  }
  if (!nfound) {    
    null = join(0,1);
    for(i=0; i<N-1; ++i) {
      null = join(null,i+1);
    }
    N += 1;
  }

  // Now the 1 element, i.e. the element of the lattice such that x ~ 1 for all x
  nfound = false;
  for(i=0; i<N; ++i) {
    found = true;
    for(j=0; j<N; ++j) {
      if (i == j) continue;
      if (get_order(i,j) != Relation::after) {
        found = false;
        break;
      }
    }
    if (found) {
      unity = i;
      nfound = true;
      break;
    }
  }
  if (!nfound) {
    unity = meet(0,1);
    for(i=0; i<N-1; ++i) {
      unity = meet(unity,i+1);
    }
    N += 1;
  }   
}

void Lattice::compute_atoms()
{
  int i,j;
  bool found;
  std::set<int>::const_iterator it;

  atoms.clear();
  for(i=0; i<N; ++i) {
    if (i == null) continue;
    if (covered(null,i)) atoms.insert(i);
  }
  
  atomic = true;
  for(i=0; i<N; ++i) {
    if (i == null) continue;
    if (atoms.count(i) > 0) continue;
    found = false;
    for(it=atoms.begin(); it!=atoms.end(); ++it) {
      j = *it;
      if (get_order(j,i) == Relation::before) {
        found = true;
        break;
      }
    }
    if (!found) {
      atomic = false;
      break;
    }
  }
}

int Lattice::meet(int x,int y) const
{
  // Find the element w in L such that w ~ x and w ~ y, while any other 
  // z in L satisfying these relations is such that z ~ w.
  int i,j,rvalue = N;
  bool max;
  std::set<int> candidates;
  std::set<int>::const_iterator it,jt;
  Relation rho;
  
  for(i=0; i<N; ++i) {
    rho = get_order(i,x);
    if (rho != Relation::before) continue;
    rho = get_order(i,y);
    if (rho != Relation::before) continue;
    candidates.insert(i); 
  }
  // Find which among the candidates is the largest
  for(it=candidates.begin(); it!=candidates.end(); ++it) {
    i = *it;
    max = true;
    for(jt=candidates.begin(); jt!=candidates.end(); ++jt) {
      j = *jt;
      if (i == j) continue;
      rho = get_order(i,j);
      if (rho == Relation::before) {
        max = false;
        break;
      }
    }
    if (max) return i;
  }
  return rvalue;
}

int Lattice::join(int x,int y) const
{
  // Find the element w in L such that x ~ w and y ~ w, while any other 
  // z in L satisfying these relations is such that w ~ z.
  int i,j,rvalue = N;
  bool max;
  std::set<int> candidates;
  std::set<int>::const_iterator it,jt;
  Relation rho;
  
  for(i=0; i<N; ++i) {
    rho = get_order(i,x);
    if (rho != Relation::after && i != x) continue;
    rho = get_order(i,y);
    if (rho != Relation::after && i != y) continue;
    candidates.insert(i); 
  }

  // Find which among the candidates is the largest
  for(it=candidates.begin(); it!=candidates.end(); ++it) {
    i = *it;
    max = true;
    for(jt=candidates.begin(); jt!=candidates.end(); ++jt) {
      j = *jt;
      if (i == j) continue;
      rho = get_order(i,j);
      if (rho == Relation::after) {
        max = false;
        break;
      }
    }
    if (max) return i;
  }
  return rvalue;
}

