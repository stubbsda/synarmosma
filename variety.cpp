/*
  Copyright 2014 Daniel Stubbs

  This file is part of Synarmosma.

  Synarmosma is free software: you can redistribute it and/or modify 
  it under the terms of the GNU General Public License as published by 
  the Free Software Foundation, either version 3 of the License, or 
  (at your option) any later version.

  Synarmosma is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with Synarmosma.  If not, see <http://www.gnu.org/licenses/>.
*/

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
Variety<kind>::Variety(unsigned int n)
{
  assert(n > 0);
  nequation = n;
  characteristic = 0;
  allocate();
}

namespace SYNARMOSMA 
{
  template<>
  Variety<unsigned int>::Variety(unsigned int n,unsigned int p)
  {
    assert(n > 0 && p > 0);
    nequation = n;
    if (!NTL::ProbPrime(p)) {
      std::cerr << "The field characteristic must be a prime!" << std::endl;
      std::exit(1);
    }
    nvariable = p;
    characteristic = p;
    allocate();
  }
}

template<class kind>
Variety<kind>::Variety(unsigned int n,unsigned int p)
{
  assert(n > 0 && p == 0);
  nequation = n;
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
  assert(nequation > 0);
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
  dependencies.clear();
  remainder.clear();
  delete[] equations;
  nequation = 0;
}

template<class kind>
void Variety<kind>::elaborate()
{
  unsigned int i,j,k;
  std::set<unsigned int> atoms;

  dependencies.clear();

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
    assert(characteristic > 0);
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
  assert(characteristic == 0);
  return 0;
}

template<class kind>
void Variety<kind>::find_partial(std::vector<unsigned int>& connected,unsigned int first,const std::vector<unsigned int>* dual_system) const
{
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
  void Variety<unsigned int>::zeta_function(unsigned int k,std::vector<unsigned int>& output) const
  {
    assert(k > 0);
    assert(characteristic > 0);
    unsigned int i,j,m,n,l,bdry,in1,f,size,nsolution;
    bool soln,found;
    Monomial<unsigned int> term;

    NTL::ZZ_p::init(NTL::to_ZZ(characteristic));

    // First we calculate the case $k = 1$
    output.clear();
    output.push_back(compute_zeros()); 

    // Now we need to calculate the case $k > 1$ 
    for(n=2; n<=k; ++n) {
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
void Variety<kind>::zeta_function(unsigned int k,std::vector<unsigned int>& output) const
{
  assert(characteristic == 0);
}


namespace SYNARMOSMA {
  template<>
  void Variety<unsigned int>::normalize(unsigned int n)
  {
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
void Variety<kind>::normalize(unsigned int n)
{
  assert(characteristic == 0);
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
void Variety<kind>::set_remainder_value(unsigned int n,kind r)
{
  assert(n < nequation);
  remainder[n] = r;
}

template<class kind>
void Variety<kind>::add_term(unsigned int n,const Monomial<kind>& t)
{
  assert(n < nequation);
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
void Variety<kind>::add_term(unsigned int n,kind alpha,const std::vector<unsigned int>& xp)
{
  assert(n < nequation);
  assert(xp.size() == nvariable);
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

