#include "variety.h"

using namespace SYNARMOSMA;

extern Random RND;

template<class kind>
Variety<kind>::Variety()
{
  nequation = 3 + RND.irandom(5);
  characteristic = 0;
  allocate();
}

template<class kind>
Variety<kind>::Variety(int p)
{
  nequation = p;
  characteristic = 0;
  allocate();
}

template<class kind>
Variety<kind>::Variety(int n,int p)
{
  nequation = n;
  nvariable = p;
  characteristic = p;
  allocate();
}

template<class kind>
Variety<kind>::~Variety()
{
  delete[] equations;
}

template<class kind>
void Variety<kind>::allocate()
{
  equations = new std::vector<Monomial<kind> >[nequation];
  initialize();
}

template<class kind>
void Variety<kind>::initialize()
{
  int i,j,k,l,alpha,beta;
  unsigned int test;
  Monomial<kind> term;
  std::set<int> atoms;
  std::pair<int,int> duo;
  bool good;
  kind in1;

  projective = false;
  for(i=0; i<nequation; ++i) {
    alpha = RND.irandom(3,8);
    atoms.clear();
    for(j=0; j<alpha; ++j) {
      beta = RND.irandom(1,5);
      term.coefficient = (characteristic > 0) ? RND.irandom(characteristic) : RND.irandom(-10,10);
      assert(beta < nvariable);
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
  int i;
  dependencies.clear();
  for(i=0; i<nequation; ++i) {
    equations[i].clear();
    remainder[i] = 0;
  }
}

template<class kind>
void Variety<kind>::elaborate()
{
  int i,j,k;
  std::set<int> atoms;

  for(i=0; i<nequation; ++i) {
    atoms.clear();
    for(j=0; j<(signed) equations[i].size(); ++j) {
      for(k=0; k<(signed) equations[i][j].exponents.size(); ++k) {
        atoms.insert(equations[i][j].exponents[k].first);
      }
    }
    dependencies.push_back(atoms);
  }
}

template<class kind>
void Variety<kind>::make_projective()
{
  if (projective) return;
  int i,j,k,in1,sum;
  Monomial<kind> term;
  bool equal;
  std::vector<int> exponents;
  std::pair<int,int> duo;

  duo.first = nvariable + 1;
  for(i=0; i<nequation; ++i) {
    for(j=0; j<(signed) equations[i].size(); ++j) {
      term = equations[i][j];
      sum = 0;
      for(k=0; k<(signed) term.exponents.size(); ++k) {
        sum += term.exponents[k].second;
      }
      exponents.push_back(sum);
    }
    // We need to know if all the values of "exponents" are equal
    equal = true;
    in1 = exponents[0];
    for(j=1; j<(signed) exponents.size(); ++j) {
      if (exponents[j] != in1) equal = false;
      if (exponents[j] > in1) in1 = exponents[j];
    }
    if (!equal) {
      for(j=0; j<(signed) equations[i].size(); ++j) {
        if (exponents[j] < in1) {
          duo.second = in1 - exponents[j];
          equations[i][j].exponents.push_back(duo);
        }
      }
    }
    exponents.clear();
    if (remainder[i] != kind(0)) {
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
  int Variety<Rational>::compute_zeros()
  {
    return 0;
  }
}

template<class kind>
int Variety<kind>::compute_zeros()
{
  int bdry,nsolution,i,j,k,l,f,in1,value,wt;
  std::vector<int> elements,vec;
  bool soln;
  Monomial<kind> term;

  assert(characteristic > 0);
  bdry = ipow(characteristic,nvariable);
  for(i=0; i<characteristic; ++i) {
    elements.push_back(i);
  }
  nsolution = 0;
  for(i=0; i<bdry; ++i) {
    in1 = i;
    for(j=0; j<nvariable; ++j) {
      //p = ipow(characteristic,nvariable-j+1);
      f = in1/characteristic;
      vec.push_back(elements[f]);
      in1 -= f*characteristic;
    }
    soln = true;
    for(j=0; j<nequation; ++j) {
      value = 0;
      for(k=0; k<(signed) equations[j].size(); ++k) {
        term = equations[j][k];
        wt = term.coefficient;
        for(l=0; l<(signed) term.exponents.size(); ++l) {
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

template<class kind>
void Variety<kind>::find_partial(bool* connected,int first,const std::vector<int>* dual_system) const
{
  int i,in1,in2;
  bool done = false;
  std::set<int>::const_iterator it,jt;
  std::set<int> current,next,handled;

  current = dependencies[first];

  do {
    for(it=current.begin(); it!=current.end(); ++it) {
      in1 = *it;
      for(i=0; i<(signed) dual_system[in1].size(); ++i) {
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
      connected[*it] = true;
    }
    handled.clear();
    if (current.empty()) done = true;
  } while(!done);
}

template<class kind>
int Variety<kind>::compute_dependencies(int* belongs) const
{
  int i,j,first = 0,family = 0;
  bool done = false;
  bool* connected = new bool[nequation];
  std::vector<int>* dual_system = new std::vector<int>[nvariable];

  for(i=0; i<nvariable; ++i) {
    // We must determine which variables this equation contains
    for(j=0; j<nequation; ++j) {
      if (dependencies[j].count(i) > 0) dual_system[i].push_back(j);
    }
  }

  for(i=0; i<nequation; ++i) {
    connected[i] = false;
  }
  do {
    done = true;
    connected[first] = true;
    find_partial(connected,first,dual_system);
    for(i=0; i<nequation; ++i) {
      if (connected[i]) {
        belongs[i] = family;
      }
      else {
        first = i;
        done = false;
      }
    }
    family += 1;
    for(i=0; i<nequation; ++i) {
      connected[i] = false;
    }
  } while(!done);

  delete[] dual_system;
  delete[] connected;
  return family;
}

namespace SYNARMOSMA {
  template<>
  void Variety<Rational>::zeta_function(int k,int* output)
  {

  }
}

template<class kind>
void Variety<kind>::zeta_function(int k,int* output)
{
  assert(k >= 1);
  assert(characteristic > 0);
  int i,j,m,n,l,bdry,in1,f,size,nsolution;
  bool soln,found;
  Monomial<kind> term;

  NTL::ZZ_p::init(NTL::to_ZZ(characteristic));


  // First we calculate the case $k = 1$
  output[0] = compute_zeros(); 
  if (k == 1) return;

  // Now we need to calculate the case $k > 1$ 
  n = 2;
  do {
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
      for(i=0; i<(signed) elements.size(); ++i) {
        if (wt == elements[i]) {
          found = true;
          break;
        }
      }
      if (!found) {
        elements.push_back(wt);
      }
      if (elements.size() == NTL::ZZ_pE::cardinality()) break;
    } while(true);

    nsolution = 0;
    for(i=0; i<bdry; ++i) {
      in1 = i;
      for(j=0; j<nvariable; ++j) {
        //p = ipow(size,nvariable-j+1);
        f = in1/size;
        vec.push_back(elements[f]);
        in1 -= f*size;
      }
      soln = true;
      for(j=0; j<nequation; ++j) {
        value = NTL::ZZ_pE::zero();
        for(m=0; m<(signed) equations[j].size(); ++m) {
          term = equations[j][m];
          wt = term.coefficient;
          for(l=0; l<(signed) term.exponents.size(); ++l) {
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
    output[n-1] = nsolution;
  } while(n <= k);
}

namespace SYNARMOSMA {
  template<>
  void Variety<int>::normalize(int n)
  {
    if (characteristic == 0) return;
    int i,in1;
    for(i=0; i<(signed) equations[n].size(); ++i) {
      in1 = equations[n][i].coefficient;
      in1 = in1 % characteristic;
      equations[n][i].coefficient = in1;
    }
  }
}

template<class kind>
void Variety<kind>::normalize(int n)
{

}

template<class kind>
void Variety<kind>::set_value(int n,kind r)
{
  remainder[n] = r;

}

template<class kind>
void Variety<kind>::add_term(int n,const Monomial<kind>& t)
{
  assert(n < nequation);
  // Have we already seen this term?
  int i,in1;
  bool found = false;
  std::vector<std::pair<unsigned int,unsigned int> > power,term = t.exponents;
  for(i=0; i<(signed) equations[n].size(); ++i) {
    power = equations[n][i].exponents;
    if (power == term) {
      in1 = i;
      found = true;
      break;
    }
  }
  if (found) {
    equations[n][in1].coefficient = equations[n][in1].coefficient + t.coefficient;
    normalize(n);
  }
  else {
    equations[n].push_back(t);
  }
}

template<class kind>
void Variety<kind>::add_term(int n,kind alpha,const int* xp)
{
  assert(n < nequation);
  int i,j,in1,sum = 0;
  Monomial<kind> term;
  bool equal;
  std::pair<int,int> duo;
  std::vector<int> exponents;

  term.coefficient = alpha;
  for(i=0; i<nvariable; ++i) {
    if (xp[i] > 0) {
      duo.first = i;
      duo.second = xp[i];
      term.exponents.push_back(duo);
      sum += xp[i];
      dependencies[n].insert(i);
    }
  }
  equations[n].push_back(term);

  if (sum > 1) linear = false;
  if (sum == 0) homogeneous = false;
  for(i=0; i<(signed) equations[n].size(); ++i) {
    term = equations[n][i];
    sum = 0;
    for(j=0; j<(signed) term.exponents.size(); ++j) {
      sum += term.exponents[j].second;
    }
    exponents.push_back(sum);
  }
  equal = true;
  in1 = exponents[0];
  for(i=1; i<(signed) exponents.size(); ++i) {
    if (exponents[i] != in1) {
      equal = false;
      break;
    }
  }
  if (!equal && projective) projective = false;
  if (equal && !projective) {
    // We need to know whether this equation was the only reason why the variety wasn't projective

  }
}

