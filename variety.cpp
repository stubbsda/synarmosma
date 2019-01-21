#include "variety.h"

using namespace SYNARMOSMA;

extern Random RND;

template<class kind>
Variety<kind>::Variety()
{

}

template<class kind>
Variety<kind>::Variety(int n)
{
  if (n < 1) throw std::invalid_argument("The number of equations in the variety must be greater than zero!");

  nequation = n;
  nvariable = n;
  characteristic = 0;
  allocate();
}

namespace SYNARMOSMA {
  template<>
  Variety<Rational>::Variety(int n,int p) 
  {
    if (n < 1) throw std::invalid_argument("The number of equations in the variety must be greater than zero!");
    if (p != 0) throw std::invalid_argument("The characteristic must be zero when the base field is rational!");
    nequation = n;
    nvariable = n;
    characteristic = 0;
    allocate();
  }
}

template<class kind>
Variety<kind>::Variety(int n,int p)
{
  if (n < 1) throw std::invalid_argument("The number of equations in the variety must be greater than zero!");
  if (p > 0) {
    if (!NTL::ProbPrime(p)) throw std::invalid_argument("The field characteristic must be zero or a prime!");
  }
  else if (p < 0) {
    throw std::invalid_argument("The field characteristic must be zero or a prime!");
  }
  nequation = n;
  nvariable = p;
  characteristic = p;
  allocate();
}

template<class kind>
Variety<kind>::Variety(const Variety<kind>& source)
{
  clear();

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

  clear();

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
  if (nequation < 1) throw std::invalid_argument("The number of equations to be allocated must be greater than zero!");

  equations = new std::vector<Monomial<kind> >[nequation];
}

namespace SYNARMOSMA {
  template<>
  Rational Variety<Rational>::generate_coefficient(int L) const
  {
    int n = RND.irandom(-L*L,L*L);
    int d = RND.irandom(-L,L);
    Rational output(n,d);
    return output;
  }

  template<>
  NTL::ZZ Variety<NTL::ZZ>::generate_coefficient(int L) const
  {
    int t = (characteristic > 0) ? RND.irandom(characteristic) : RND.irandom(-L,L);
    NTL::ZZ output = NTL::to_ZZ(t);
    return output;
  }
}

template<class kind>
kind Variety<kind>::generate_coefficient(int L) const
{
  kind output = (characteristic > 0) ? RND.irandom(characteristic) : RND.irandom(-L,L);
  return output;
}

template<class kind>
void Variety<kind>::random_variety(unsigned int mdegree)
{
  if (mdegree < 1) throw std::invalid_argument("The maximum degree of a randomly generated variety must be greater than zero!");
  unsigned int i,j,k,l,nterm,alpha,beta,test;
  Monomial<kind> term;
  std::set<unsigned int> atoms;
  std::pair<unsigned int,unsigned int> duo;
  bool good;
  kind rho;

  projective = false;
  for(i=0; i<nequation; ++i) {
    nterm = RND.irandom(3,8);
    if (nterm > mdegree*nvariable) nterm = mdegree*nvariable;
    atoms.clear();
    for(j=0; j<nterm; ++j) {
      beta = RND.irandom(1,1+nvariable);
      term.coefficient = generate_coefficient(10);
      tpower = 0;
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
        alpha = RND.irandom(1,1+mdegree);
        duo.second = ((tpower + alpha) <= mdegree) ? alpha : 0;
        if (duo.second > 0) {
          term.exponents.push_back(duo);
          atoms.insert(duo.first);
          tpower += duo.second;
        }
        else {
          break;
        }
      }
      equations[i].push_back(term);
      term.exponents.clear();
    }
    rho = generate_coefficient(10);
    remainder.push_back(rho);
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
        count += write_ZZ(s,equations[i][j].coefficient.get_numerator());
        count += write_ZZ(s,equations[i][j].coefficient.get_denominator());
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
      count += write_ZZ(s,remainder[i].get_numerator());
      count += write_ZZ(s,remainder[i].get_denominator());
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
void Variety<kind>::make_projective()
{
  if (projective) return;
  unsigned int i,j,k,in1,sum;
  Monomial<kind> term;
  bool equal;
  std::vector<unsigned int> exponents;
  std::pair<unsigned int,unsigned int> duo;

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
      remainder[i] = zero;
    }
  }
  nvariable += 1;
  projective = true;
  homogeneous = true;
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

template<class kind>
void Variety<kind>::compute_properties()
{
  unsigned int i,j,k,sum;
  std::set<unsigned int> tpower;

  linear = true;
  homogeneous = true;
  projective = true;
  for(i=0; i<nequation; ++i) {
    if (remainder[i] != zero) {
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
bool Variety<kind>::add_term(int n,const Monomial<kind>& t)
{
  if (n < 0 || n >= (signed) nequation) throw std::invalid_argument("Illegal equation number in Variety::add_term!");

  // Have we already seen this term?
  unsigned int i,in1;
  bool found = false,output = true;
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
    output = false;
  }
  else {
    equations[n].push_back(t);
  }
  compute_properties();

  return output;
}

template<class kind>
bool Variety<kind>::add_term(int n,kind alpha,const std::vector<unsigned int>& powers)
{
  if (powers.size() != nvariable) throw std::invalid_argument("The length of the vector of exponents must be equal to the number of variables in the variety!");

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
  return add_term(n,term);
}

