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

#ifndef _globalh
#define _globalh

// STL headers
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cassert>
#include <vector>
#include <cmath>
#include <set>
#include <algorithm>
#include <ctime>
#include <string>
#include <limits>
#include <sys/time.h>
// NTL headers for number theory functionality...
#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/GF2.h>
// The Boost library headers...
#include <boost/tokenizer.hpp>
#include <boost/date_time.hpp>
#include <boost/timer/timer.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/unordered_map.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/poisson_distribution.hpp>
#include <boost/random/bernoulli_distribution.hpp>
#include <boost/math/distributions/beta.hpp>
// Finally the LAPACK declarations
extern "C" {
  void dgemm_(char*,char*,int*,int*,int*,double*,double*,int*,double*,int*,double*,double*,int*);
  void dsyev_(char*,char*,int*,double*,int*,double*,double*,int*,int*);
  void dgetri_(int*,double*,int*,int*,double*,int*,int*);
  void dgetrf_(int*,int*,double*,int*,int*,int*);
}

namespace SYNARMOSMA {
  typedef boost::mt19937 base_generator_type;
  typedef boost::unordered_map<std::set<int>,int> hash_map;
  typedef boost::unordered_map<std::string,int> string_hash;
  typedef boost::adjacency_list<boost::vecS,boost::vecS,boost::bidirectionalS> network;

#ifdef __LP64__
  typedef unsigned long UINT64;
#else
  typedef unsigned long long UINT64;
#endif

  template<unsigned int N>
  struct power_of_two {
    static unsigned int const value = 2*power_of_two<N-1>::value;
  };

  template<>
  struct power_of_two<0> {
    static unsigned int const value = 1;
  };

  template<class kind>
  class Monomial {
   public:
    kind coefficient;
    std::vector<std::pair<unsigned int,unsigned int> > exponents; 
  };

  UINT64 ipow(int,int);
  UINT64 factorial(int);
  int combinations(const std::set<int>&,int,std::vector<int>&);
  void get_neighbours(int,const std::vector<int>*,int,std::vector<int>&);
  void induced_orientation(int,const std::vector<int>&,int,const hash_map&,int*);
  void factorize(long,std::vector<std::pair<long,int> >&);
  int parity(const std::vector<int>&,const std::vector<int>&);
  void invert(const double*,double*,int);
  double determinant(const double*,int);
  void convert(unsigned char*,int);
  void convert(unsigned char*,float);
  void compute_smith_normal_form(std::vector<std::pair<int,int> >*,int,int);
  void vertex_difference(int,int,std::vector<double>&);
  double binomial(int,int);
  double arithmetic_mean(const std::vector<double>&);
  void split(const std::string&,char,std::vector<std::string>&);
  void split(const std::vector<int>&,std::vector<int>&);
  double dmap(double,double);
  double norm(const double*,int);
  double norm(const std::vector<double>&);
  void cross_product(const std::vector<double>&,const std::vector<double>&,std::vector<double>&);
  bool tuple_predicate(const boost::tuple<int,int,double>&,const boost::tuple<int,int,double>&);
  bool pair_predicate_dbl(const std::pair<int,double>&,const std::pair<int,double>&);
  bool pair_predicate_int(const std::pair<int,int>&,const std::pair<int,int>&);
  int element(const std::set<int>&);
  bool next_combination(std::vector<int>&,int);
  void complement(const std::set<int>&,std::set<int>&,int,int);
  int coincidence(const std::set<int>&,const std::set<int>&);
  bool double_equality(double,double);
  void trim(std::string&);
  void set_wcomponent_values(int,bool);

  inline void RGB_intensity(double rho,unsigned char* output)  
  {
    // This assumes that $\rho \in [0,\pi/2]$
    double x = std::sin(rho);
    double y = std::cos(rho);
    double z = std::sin(2.0*rho);
    // Red
    output[0] = (unsigned char) (255.0*x*x);
    // Green
    output[1] = (unsigned char) (255.0*z);
    // Blue
    output[2] = (unsigned char) (255.0*y*y);
  }

  inline std::string make_key(int x)
  {
    std::stringstream s;
    s << x;
    return s.str();
  }

  inline std::string make_key(int x,int y)
  {
    assert(x != y);
    std::stringstream s;
    if (x < y) {
      s << x << ":" << y;
    }
    else {
      s << y << ":" << x;
    }
    return s.str();
  }

  inline std::string make_key(const std::vector<int>& v) 
  {
    assert(!v.empty());
    unsigned int i,n = v.size();
    std::stringstream s;

    for(i=0; i<n-1; ++i) {
      s << v[i] << ":";
    }
    s << v[n-1];
    return s.str();
  }

  inline std::string make_key(const std::set<int>& S)
  {
    assert(!S.empty());
    unsigned int i,n = S.size();
    std::stringstream s;
    std::set<int>::const_iterator it;
    i = 0;
    for(it=S.begin(); it!=S.end(); ++it) {
      if (i < (n-1)) {
        s << *it << ":";
      }
      else {
        s << *it;
      }
      ++i;
    }
    return s.str();
  }

  class Random {
   private:
    unsigned int s;
    base_generator_type BGT;

    boost::math::beta_distribution<>* root_beta;

    boost::bernoulli_distribution<>* brn;
    boost::variate_generator<base_generator_type&,boost::bernoulli_distribution<> >* vbrn;

    boost::poisson_distribution<>* poisson;
    boost::variate_generator<base_generator_type&,boost::poisson_distribution<> >* vpoisson;

    boost::uniform_real<>* uniform;
    boost::variate_generator<base_generator_type&,boost::uniform_real<> >* VRG;

    boost::normal_distribution<>* gaussian;
    boost::variate_generator<base_generator_type&,boost::normal_distribution<> >* NRG;

    bool brn_allocated;
    bool beta_allocated;
   public:
    Random();
    ~Random();
    inline void set_seed(unsigned int x) {s = x; BGT.seed(s);};
    inline void increment_seed() {s++; BGT.seed(s);};
    inline void decrement_seed() {s--; BGT.seed(s);};
    inline unsigned int get_seed() const {return s;};
    void initialize_beta(double,double);
    void initialize_bernoulli(double);
    void initialize_poisson(double);
    double drandom();
    double drandom(double,double);
    double beta_variate();
    bool bernoulli_variate();
    bool poisson_variate();
    double nrandom();
    double nrandom(double);
    double nrandom(double,double);
    int irandom(int);
    int irandom(int,int);
    int irandom(const std::set<int>&);
    int irandom(const std::set<int>&,const std::set<int>&);
    int irandom(const std::vector<int>&);
    int irandom(int,const std::vector<int>&);
    void shuffle(std::vector<int>&,int);
  };
}
#endif
