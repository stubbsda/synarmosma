#include "variety.h"

using namespace SYNARMOSMA;

extern Random RND;

template<class kind>
Variety<kind>::Variety()
{
  nequation = 3 + RND.irandom(5);
  nvariable = 2 + RND.irandom(3);
  characteristic = 0;
  allocate();
}

template<class kind>
Variety<kind>::Variety(int n)
{
#ifdef DEBUG
  assert(n > 0);
#endif
  nequation = n;
  nvariable = n;
  characteristic = 0;
  allocate();
}

namespace SYNARMOSMA 
{
  template<>
  Variety<unsigned int>::Variety(int n,int p)
  {
#ifdef DEBUG
    assert(n > 0 && p > 0);
#endif
    nequation = n;
    nvariable = n;
    if (!NTL::ProbPrime(p)) throw std::invalid_argument("Field characteristic must be prime!");
    nvariable = p;
    characteristic = p;
    allocate();
  }
}

template<class kind>
Variety<kind>::Variety(int n,int p)
{
#ifdef DEBUG
  assert(n > 0 && p == 0);
#endif
  nequation = n;
  nvariable = n;
  nvariable = p;
  characteristic = p;
  allocate();
}

template<class kind>
Variety<kind>::Variety(const Variety<kind>& source)
{
  unsigned int i;
  remainder = source.remainder;
  nequation = source.nequation;
  nvariable = source.nvariable;
  characteristic = source.characteristic;
  linear = source.linear;
  homogeneous = source.homogeneous;
  projective = source.projective;
  equations = new std::vector<Monomial<kind> >[nequation];
  for(i=0; i<nequation; ++i) {
    equations[i] = source.equations[i];
  }
}

template<class kind>
Variety<kind>& Variety<kind>::operator =(const Variety<kind>& source)
{
  if (this == &source) return *this;

  if (nequation > 0) delete[] equations;
  unsigned int i;
  remainder = source.remainder;
  nequation = source.nequation;
  nvariable = source.nvariable;
  characteristic = source.characteristic;
  linear = source.linear;
  homogeneous = source.homogeneous;
  projective = source.projective;
  equations = new std::vector<Monomial<kind> >[nequation];
  for(i=0; i<nequation; ++i) {
    equations[i] = source.equations[i];
  }

  return *this;
}

template<class kind>
Variety<kind>::~Variety()
{
  if (nequation > 0) delete[] equations;
}

template<class kind>
void Variety<kind>::allocate()
{
#ifdef DEBUG
  assert(nequation > 0);
#endif
  equations = new std::vector<Monomial<kind> >[nequation];
  initialize();
}

template<class kind>
void Variety<kind>::initialize()
{
  unsigned int i,j,k,l,alpha,beta,test;
  Monomial<kind> term;
  std::set<unsigned int> atoms;
  std::pair<unsigned int,unsigned int> duo;
  bool good;
  kind in1;

  projective = false;
  for(i=0; i<nequation; ++i) {
    alpha = RND.irandom(3,8);
    atoms.clear();
    for(j=0; j<alpha; ++j) {
      beta = RND.irandom(1,1+nvariable);
      term.coefficient = (characteristic > 0) ? RND.irandom(characteristic) : RND.irandom(-10,10);
#ifdef DEBUG
      assert(beta <= nvariable);
#endif
      for(k=0; k<beta; ++k) {
        // What of the case where duo.first is the same number for different iterations
        // of this loop over k? We need to ensure that this never happens...
        do {
          test = RND.irandom(nvariable);
          good = true;
          for(l=0; l<k; ++l) {
            if (term.exponents[l].first == test) {
              good = false;
              break;
            }
          }
          if (good) break;
        } while(true);
        duo.first = test;
        duo.second = RND.irandom(1,6);
        term.exponents.push_back(duo);
        atoms.insert(duo.first);
      }
      equations[i].push_back(term);
      term.exponents.clear();
      in1 = (characteristic > 0) ? RND.irandom(characteristic) : RND.irandom(-10,10);
      remainder.push_back(in1);
    }
    dependencies.push_back(atoms);
  }
}

template<class kind>
void Variety<kind>::clear()
{
  // This method gets called by the Functional_Equation::reduce as the final step
  // in building up a variety over GF(p) from reducing a functional equation over this
  // same finite field. The method's job is to build up the contents of the vector of
  // sets containing the variable dependency information for each equation in this
  // algebraic variety over GF(p).
  dependencies.clear();
  remainder.clear();
  if (nequation > 0) delete[] equations;
  nequation = 0;
  nvariable = 0;
  characteristic = 0;
  linear = false;
  homogeneous = false;
  projective = false;
}

namespace SYNARMOSMA {
  template<>
  int Variety<Rational>::write_equations(std::ofstream& s) const
  {
    unsigned int i,j,k,l,n,m;
    int count = 0;

    for(i=0; i<nequation; ++i) {
      n = equations[i].size();
      s.write((char*)(&n),sizeof(int)); count += sizeof(int);
      for(j=0; j<n; ++j) {
        count += write_ZZ(s,equations[i][j].coefficient.numerator());
        count += write_ZZ(s,equations[i][j].coefficient.denominator());
        m = equations[i][j].exponents.size();
        s.write((char*)(&m),sizeof(int)); count += sizeof(int);
        for(k=0; k<m; ++k) {
          l = equations[i][j].exponents[k].first;
          s.write((char*)(&l),sizeof(int)); count += sizeof(int);
          l = equations[i][j].exponents[k].second;
          s.write((char*)(&l),sizeof(int)); count += sizeof(int);
        }
      }
    }
    for(i=0; i<nequation; ++i) {
      count += write_ZZ(s,remainder[i].numerator());
      count += write_ZZ(s,remainder[i].denominator());
    }
    return count;
  }

  template<>
  int Variety<NTL::ZZ>::write_equations(std::ofstream& s) const
  {
    unsigned int i,j,k,l,n,m;
    int count = 0;

    for(i=0; i<nequation; ++i) {
      n = equations[i].size();
      s.write((char*)(&n),sizeof(int)); count += sizeof(int);
      for(j=0; j<n; ++j) {
        write_ZZ(s,equations[i][j].coefficient);
        m = equations[i][j].exponents.size();
        s.write((char*)(&m),sizeof(int)); count += sizeof(int);
        for(k=0; k<m; ++k) {
          l = equations[i][j].exponents[k].first;
          s.write((char*)(&l),sizeof(int)); count += sizeof(int);
          l = equations[i][j].exponents[k].second;
          s.write((char*)(&l),sizeof(int)); count += sizeof(int);
        }
      }
    }
    for(i=0; i<nequation; ++i) {
      count += write_ZZ(s,remainder[i]);
    }
    return count;
  }
}

template<class kind>
int Variety<kind>::write_equations(std::ofstream& s) const
{
  unsigned int i,j,k,l,n,m;
  int count = 0;
  kind x;

  for(i=0; i<nequation; ++i) {
    n = equations[i].size();
    s.write((char*)(&n),sizeof(int)); count += sizeof(int);
    for(j=0; j<n; ++j) {
      x = equations[i][j].coefficient;
      s.write((char*)(&x),sizeof(kind)); count += sizeof(kind);
      m = equations[i][j].exponents.size();
      s.write((char*)(&m),sizeof(int)); count += sizeof(int);
      for(k=0; k<m; ++k) {
        l = equations[i][j].exponents[k].first;
        s.write((char*)(&l),sizeof(int)); count += sizeof(int);
        l = equations[i][j].exponents[k].second;
        s.write((char*)(&l),sizeof(int)); count += sizeof(int);
      }
    }
  }
  for(i=0; i<nequation; ++i) {
    x = remainder[i];
    s.write((char*)(&x),sizeof(kind)); count += sizeof(kind);
  }
  return count;
}

template<class kind>
int Variety<kind>::serialize(std::ofstream& s) const
{
  int count = 0;

  s.write((char*)(&nequation),sizeof(int)); count += sizeof(int);
  s.write((char*)(&nvariable),sizeof(int)); count += sizeof(int);
  s.write((char*)(&characteristic),sizeof(int)); count += sizeof(int);
  s.write((char*)(&linear),sizeof(bool)); count += sizeof(bool);
  s.write((char*)(&homogeneous),sizeof(bool)); count += sizeof(bool);
  s.write((char*)(&projective),sizeof(bool)); count += sizeof(bool);
  count += write_equations(s);

  return count;
}

namespace SYNARMOSMA {
  template<>
  int Variety<Rational>::read_equations(std::ifstream& s)
  {
    unsigned int i,j,k,l1,l2,n,m;
    int count = 0;
    NTL::ZZ q1,q2;
    Monomial<Rational> t;

    for(i=0; i<nequation; ++i) {
      s.read((char*)(&n),sizeof(int)); count += sizeof(int);
      for(j=0; j<n; ++j) {
        count += read_ZZ(s,q1);
        count += read_ZZ(s,q2);
        t.coefficient = Rational(q1,q2);
        s.read((char*)(&m),sizeof(int)); count += sizeof(int);
        for(k=0; k<m; ++k) {
          s.read((char*)(&l1),sizeof(int)); count += sizeof(int);
          s.read((char*)(&l2),sizeof(int)); count += sizeof(int);
          t.exponents.push_back(std::pair<unsigned int,unsigned int>(l1,l2));
        }
        equations[i].push_back(t);
        t.exponents.clear();
      }
    }
    for(i=0; i<nequation; ++i) {
      count += read_ZZ(s,q1);
      count += read_ZZ(s,q2);
      remainder.push_back(Rational(q1,q2));
    }
    return count;
  }

  template<>
  int Variety<NTL::ZZ>::read_equations(std::ifstream& s)
  {
    unsigned int i,j,k,l1,l2,n,m;
    int count = 0;
    NTL::ZZ q;
    Monomial<NTL::ZZ> t;

    for(i=0; i<nequation; ++i) {
      s.read((char*)(&n),sizeof(int)); count += sizeof(int);
      for(j=0; j<n; ++j) {
        count += read_ZZ(s,q);
        t.coefficient = q;
        s.read((char*)(&m),sizeof(int)); count += sizeof(int);
        for(k=0; k<m; ++k) {
          s.read((char*)(&l1),sizeof(int)); count += sizeof(int);
          s.read((char*)(&l2),sizeof(int)); count += sizeof(int);
          t.exponents.push_back(std::pair<unsigned int,unsigned int>(l1,l2));
        }
        equations[i].push_back(t);
        t.exponents.clear();
      }
    }
    for(i=0; i<nequation; ++i) {
      count += read_ZZ(s,q);
      remainder.push_back(q);
    }
    return count;
  }
}

template<class kind>
int Variety<kind>::read_equations(std::ifstream& s)
{
  unsigned int i,j,k,l1,l2,n,m;
  int count = 0;
  kind x;
  Monomial<kind> t;

  for(i=0; i<nequation; ++i) {
    s.read((char*)(&n),sizeof(int)); count += sizeof(int);
    for(j=0; j<n; ++j) {
      s.read((char*)(&x),sizeof(kind)); count += sizeof(kind);
      t.coefficient = x;
      s.read((char*)(&m),sizeof(int)); count += sizeof(int);
      for(k=0; k<m; ++k) {
        s.read((char*)(&l1),sizeof(int)); count += sizeof(int);
        s.read((char*)(&l2),sizeof(int)); count += sizeof(int);
        t.exponents.push_back(std::pair<unsigned int,unsigned int>(l1,l2));
      }
      equations[i].push_back(t);
      t.exponents.clear();
    }
  }
  for(i=0; i<nequation; ++i) {
    s.read((char*)(&x),sizeof(kind)); count += sizeof(kind);
    remainder.push_back(x);
  }
  return count;
}

template<class kind>
int Variety<kind>::deserialize(std::ifstream& s) 
{
  int count = 0;

  clear();

  s.read((char*)(&nequation),sizeof(int)); count += sizeof(int);
  s.read((char*)(&nvariable),sizeof(int)); count += sizeof(int);
  s.read((char*)(&characteristic),sizeof(int)); count += sizeof(int); 
  s.read((char*)(&linear),sizeof(bool)); count += sizeof(bool);
  s.read((char*)(&homogeneous),sizeof(bool)); count += sizeof(bool);
  s.read((char*)(&projective),sizeof(bool)); count += sizeof(bool);
 
  if (nequation > 0) {
    equations = new std::vector<Monomial<kind> >[nequation];
    count += read_equations(s);
  }

  return count;
}

template<class kind>
void Variety<kind>::elaborate()
{
  unsigned int i,j,k;
  std::set<unsigned int> atoms;

  dependencies.clear();

  for(i=0; i<nequation; ++i) {
    atoms.clear();
    for(j=0; j<equations[i].size(); ++j) {
      for(k=0; k<equations[i][j].exponents.size(); ++k) {
        atoms.insert(equations[i][j].exponents[k].first);
      }
    }
    dependencies.push_back(atoms);
  }
}

namespace SYNARMOSMA {
  template<>
  void Variety<NTL::ZZ>::make_projective()
  {
    if (projective) return;
    unsigned int i,j,k,in1,sum;
    Monomial<NTL::ZZ> term;
    bool equal;
    std::vector<unsigned int> exponents;
    std::pair<unsigned int,unsigned int> duo;
    const NTL::ZZ zero = NTL::to_ZZ(0);

    duo.first = nvariable + 1;
    for(i=0; i<nequation; ++i) {
      for(j=0; j<equations[i].size(); ++j) {
        term = equations[i][j];
        sum = 0;
        for(k=0; k<term.exponents.size(); ++k) {
          sum += term.exponents[k].second;
        }
        exponents.push_back(sum);
      }
      // We need to know if all the values of "exponents" are equal
      equal = true;
      in1 = exponents[0];
      for(j=1; j<exponents.size(); ++j) {
        if (exponents[j] != in1) equal = false;
        if (exponents[j] > in1) in1 = exponents[j];
      }
      if (!equal) {
        for(j=0; j<equations[i].size(); ++j) {
          if (exponents[j] < in1) {
            duo.second = in1 - exponents[j];
            equations[i][j].exponents.push_back(duo);
          }
        }
      }
      exponents.clear();
      if (remainder[i] != zero) {
        term.exponents.clear();
        term.coefficient = remainder[i];
        duo.second = in1;
        term.exponents.push_back(duo);
        equations[i].push_back(term);
        remainder[i] = 0;
      }
    }
    nvariable += 1;
    projective = true;
    homogeneous = true;
  }
}

template<class kind>
void Variety<kind>::make_projective()
{
  if (projective) return;
  unsigned int i,j,k,in1,sum;
  Monomial<kind> term;
  bool equal;
  std::vector<unsigned int> exponents;
  std::pair<unsigned int,unsigned int> duo;
  const kind zero = 0;

  duo.first = nvariable + 1;
  for(i=0; i<nequation; ++i) {
    for(j=0; j<equations[i].size(); ++j) {
      term = equations[i][j];
      sum = 0;
      for(k=0; k<term.exponents.size(); ++k) {
        sum += term.exponents[k].second;
      }
      exponents.push_back(sum);
    }
    // We need to know if all the values of "exponents" are equal
    equal = true;
    in1 = exponents[0];
    for(j=1; j<exponents.size(); ++j) {
      if (exponents[j] != in1) equal = false;
      if (exponents[j] > in1) in1 = exponents[j];
    }
    if (!equal) {
      for(j=0; j<equations[i].size(); ++j) {
        if (exponents[j] < in1) {
          duo.second = in1 - exponents[j];
          equations[i][j].exponents.push_back(duo);
        }
      }
    }
    exponents.clear();
    if (remainder[i] != zero) {
      term.exponents.clear();
      term.coefficient = remainder[i];
      duo.second = in1;
      term.exponents.push_back(duo);
      equations[i].push_back(term);
      remainder[i] = 0;
    }
  }
  nvariable += 1;
  projective = true;
  homogeneous = true;
}

namespace SYNARMOSMA {
  template<>
  int Variety<unsigned int>::compute_zeros() const
  {
#ifdef DEBUG
    assert(characteristic > 0);
#endif
    unsigned int bdry,nsolution,i,j,k,l,f,in1;
    std::vector<unsigned int> elements,vec;
    bool soln;
    unsigned int wt,value;
    Monomial<unsigned int> term;
  
    bdry = ipow(characteristic,nvariable);
    for(i=0; i<characteristic; ++i) {
      elements.push_back(i);
    }
    nsolution = 0;
    for(i=0; i<bdry; ++i) {
      in1 = i;
      for(j=0; j<nvariable; ++j) {      
        f = in1/characteristic;
        vec.push_back(elements[f]);
        in1 -= f*characteristic;
      }
      soln = true;
      for(j=0; j<nequation; ++j) {
        value = 0;
        for(k=0; k<equations[j].size(); ++k) {
          term = equations[j][k];
          wt = term.coefficient;
          for(l=0; l<term.exponents.size(); ++l) {
            wt *= ipow(vec[term.exponents[l].first],term.exponents[l].second);
          }
          value += wt;
        }
        value = value % characteristic;
        if (value != 0) {
          soln = false;
          break;
        }
      }
      if (soln) nsolution += 1;
      vec.clear();
    }
    return nsolution;
  }
}

template<class kind>
int Variety<kind>::compute_zeros() const
{
#ifdef DEBUG
  assert(characteristic == 0);
#endif
  return 0;
}

template<class kind>
void Variety<kind>::find_partial(std::vector<unsigned int>& connected,int first,const std::vector<unsigned int>* dual_system) const
{
  assert(first >= 0);
  unsigned int i,in1,in2;
  bool done = false;
  std::set<unsigned int>::const_iterator it,jt;
  std::set<unsigned int> current,next,handled;

  current = dependencies[first];

  do {
    for(it=current.begin(); it!=current.end(); ++it) {
      in1 = *it;
      for(i=0; i<dual_system[in1].size(); ++i) {
        in2 = dual_system[in1][i];
        if (!connected[in2]) {
          handled.insert(in2);
          for(jt=dependencies[in2].begin(); jt!=dependencies[in2].end(); ++jt) {
            next.insert(*jt);
          }
        }
      }
    }
    current = next;
    next.clear();
    for(it=handled.begin(); it!=handled.end(); ++it) {
      connected[*it] = 1;
    }
    handled.clear();
    if (current.empty()) done = true;
  } while(!done);
}

template<class kind>
int Variety<kind>::compute_dependencies(std::vector<unsigned int>& component) const
{
  unsigned int i,j,first = 0,family = 0;
  bool done = false;
  std::vector<unsigned int> connected;
  std::vector<unsigned int>* dual_system = new std::vector<unsigned int>[nvariable];

  component.clear();
  for(i=0; i<nequation; ++i) {
    component.push_back(0);
    connected.push_back(0);
  }
  for(i=0; i<nvariable; ++i) {
    // We must determine which variables this equation contains
    for(j=0; j<nequation; ++j) {
      if (dependencies[j].count(i) > 0) dual_system[i].push_back(j);
    }
  }

  do {
    done = true;
    connected[first] = 1;
    find_partial(connected,first,dual_system);
    for(i=0; i<nequation; ++i) {
      if (connected[i] == 1) {
        component[i] = family;
      }
      else {
        first = i;
        done = false;
      }
    }
    family += 1;
    for(i=0; i<nequation; ++i) {
      connected[i] = 0;
    }
  } while(!done);

  delete[] dual_system;
  return family;
}

namespace SYNARMOSMA {
  template<>
  void Variety<unsigned int>::zeta_function(int k,std::vector<unsigned int>& output) const
  {
#ifdef DEBUG
    assert(k > 0);
    assert(characteristic > 0);
#endif
    unsigned int i,j,m,n,l,bdry,in1,f,size,nsolution,ku = k;
    bool soln,found;
    Monomial<unsigned int> term;

    NTL::ZZ_p::init(NTL::to_ZZ(characteristic));

    // First we calculate the case $k = 1$
    output.clear();
    output.push_back(compute_zeros()); 

    // Now we need to calculate the case $k > 1$ 
    for(n=2; n<=ku; ++n) {
      size = ipow(characteristic,n);

      NTL::ZZ_pX P;
      NTL::BuildIrred(P,n);

      NTL::ZZ_pE::init(P);
      NTL::ZZ_pE value,wt;
      std::vector<NTL::ZZ_pE> elements,vec;

      bdry = ipow(size,nvariable);
      do {
        wt = NTL::random_ZZ_pE();
        found = false;
        for(i=0; i<elements.size(); ++i) {
          if (wt == elements[i]) {
            found = true;
            break;
          }
        }
        if (!found) elements.push_back(wt);      
        if (elements.size() == NTL::ZZ_pE::cardinality()) break;
      } while(true);

      nsolution = 0;
      for(i=0; i<bdry; ++i) {
        in1 = i;
        for(j=0; j<nvariable; ++j) {
          f = in1/size;
          vec.push_back(elements[f]);
          in1 -= f*size;
        }
        soln = true;
        for(j=0; j<nequation; ++j) {
          value = NTL::ZZ_pE::zero();
          for(m=0; m<equations[j].size(); ++m) {
            term = equations[j][m];
            wt = term.coefficient;
            for(l=0; l<term.exponents.size(); ++l) {
              wt *= power(vec[term.exponents[l].first],term.exponents[l].second);
            }
            value += wt;
          }
          if (value != NTL::ZZ_pE::zero()) {
            soln = false;
            break;
          }
        }
        if (soln) nsolution += 1;
        vec.clear();
      }
      output.push_back(nsolution);
    }
  }
}

template<class kind>
void Variety<kind>::zeta_function(int k,std::vector<unsigned int>& output) const
{
#ifdef DEBUG
  assert(k > 0);
  assert(characteristic == 0);
#endif
}


namespace SYNARMOSMA {
  template<>
  void Variety<unsigned int>::normalize(int n)
  {
    assert(n >= 0);
    if (characteristic == 0) return;
    unsigned int i,in1;
 
    for(i=0; i<equations[n].size(); ++i) {
      in1 = equations[n][i].coefficient;
      in1 = in1 % characteristic;
      equations[n][i].coefficient = in1;
    }
  }
}

template<class kind>
void Variety<kind>::normalize(int n)
{
#ifdef DEBUG
  assert(n >= 0);
  assert(characteristic == 0);
#endif
  return;
}

template<class kind>
void Variety<kind>::property_check()
{
  unsigned int i,j,k,sum;
  std::set<unsigned int> tpower;

  for(i=0; i<nequation; ++i) {
    normalize(i);
  }

  linear = true;
  homogeneous = true;
  projective = true;
  for(i=0; i<nequation; ++i) {
    if (remainder[i] != 0) {
      homogeneous = false;
      projective = false;
    }
    for(j=0; j<equations[i].size(); ++j) {
      sum = 0;
      for(k=0; k<equations[i][j].exponents.size(); ++k) {
        sum += equations[i][j].exponents[k].second;
      }
      if (sum > 1) linear = false;
      tpower.insert(sum);
    }
    if (tpower.size() > 1) projective = false;
    tpower.clear();
    if (!projective && !linear && !homogeneous) break;
  }
  elaborate();
}

template<class kind>
void Variety<kind>::set_remainder_value(int n,kind r)
{
#ifdef DEBUG
  assert(n >= 0 && n < (signed) nequation);
#endif
  remainder[n] = r;
}

template<class kind>
void Variety<kind>::add_term(int n,const Monomial<kind>& t)
{
#ifdef DEBUG
  assert(n >= 0 && n < (signed) nequation);
#endif
  // Have we already seen this term?
  unsigned int i,in1;
  bool found = false;
  std::vector<std::pair<unsigned int,unsigned int> > power,term = t.exponents;
  for(i=0; i<equations[n].size(); ++i) {
    power = equations[n][i].exponents;
    if (power == term) {
      in1 = i;
      found = true;
      break;
    }
  }
  if (found) {
    equations[n][in1].coefficient = equations[n][in1].coefficient + t.coefficient;
  }
  else {
    equations[n].push_back(t);
  }
  property_check();
}

template<class kind>
void Variety<kind>::add_term(int n,kind alpha,const std::vector<unsigned int>& xp)
{
#ifdef DEBUG
  assert(n >= 0 && n < (signed) nequation);
  assert(xp.size() == nvariable);
#endif
  unsigned int i;
  Monomial<kind> term;
  std::pair<unsigned int,unsigned int> duo;

  term.coefficient = alpha;
  for(i=0; i<nvariable; ++i) {
    if (xp[i] > 0) {
      duo.first = i;
      duo.second = xp[i];
      term.exponents.push_back(duo);
    }
  }
  add_term(n,term);
}

