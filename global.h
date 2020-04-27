#ifndef _globalh
#define _globalh

// STL headers
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <set>
#include <complex>
#include <algorithm>
#include <ctime>
#include <string>
#include <limits>
#include <tuple>
#include <queue>
#include <sys/time.h>
// NTL headers for number theory functionality...
#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/GF2.h>
#include <NTL/ZZ_pEX.h>
#include <NTL/ZZ_pXFactoring.h>
// The Boost library headers...
#include <boost/tokenizer.hpp>
#include <boost/date_time.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/unordered_map.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/boyer_myrvold_planar_test.hpp>
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
  void dgesv_(int*,int*,double*,int*,int*,double*,int*,int*);
  void zgesv_(int*,int*,std::complex<double>*,int*,int*,std::complex<double>*,int*,int*);
  void dgetri_(int*,double*,int*,int*,double*,int*,int*);
}

namespace SYNARMOSMA {
  using hash_map = boost::unordered_map<std::set<int>,int>;
  using pair_index = boost::unordered_map<std::pair<int,int>,int>;

#ifdef __LP64__
  using INT64 = signed long;
  using UINT64 = unsigned long;
#else
  using INT64 = signed long long;
  using UINT64 = unsigned long long;
#endif

#ifdef DISCRETE
  const double energy_quantum = 1e-7;
  const double space_quantum = 1e-8;
#endif

  template<class kind>
  /// A template class representing a monomial, i.e. an expression of the form \f$a \prod_{k=1}^m x_{j_k}^{n_k}\f$, with \f$j_k\f$ and \f$n_k\f$ non-negative integers.
  class Monomial {
   public:
    /// The coefficient \f$a\f$, drawn from the base class of this 
    /// template class.
    kind coefficient;
    /// The variable term, \f$\prod_{k=1}^m x_{j_k}^{n_k}\f$, here stored 
    /// as a vector of pairs of unsigned integers, the first storing the 
    /// index of the variable, thus \f$j_k\f$, while the second element of 
    /// the pair is the exponent, \f$n_k\f$.
    std::vector<std::pair<unsigned int,unsigned int> > exponents; 
  };

  enum class Relation
  {
      before,
      after,
      disparate
  };

  enum class Boolean
  {
      AND,
      OR,
      XOR
  };

  UINT64 ipow(int,int);
  UINT64 factorial(int);
  UINT64 binomial(int,int);
  int combinations(const std::set<int>&,int,std::vector<int>&);
  void permute(std::vector<std::vector<int> >&,std::vector<int>&,const std::vector<int>&,int);
  void assemble_neighbours(int,const std::vector<int>*,int,std::vector<int>&);
  void induced_orientation(int,const std::vector<int>&,int,const hash_map&,int*);
  void factorize(long,std::vector<std::pair<long,int> >&);
  int parity(const std::vector<int>&,const std::vector<int>&);
  int parity(const std::vector<int>&);
  void compute_smith_normal_form(std::vector<std::pair<int,int> >*,int,int);
  void vertex_difference(int,int,std::vector<double>&);
  double arithmetic_mean(const std::vector<double>&);
  void split(const std::string&,char,std::vector<std::string>&);
  void split(const std::vector<int>&,std::vector<int>&);
  double dmap(double,double);
  int write_ZZ(std::ofstream&,const NTL::ZZ&);
  int read_ZZ(std::ifstream&,NTL::ZZ&);
  double norm(const double*,int);
  double norm(const std::vector<double>&);
  void cross_product(const std::vector<double>&,const std::vector<double>&,std::vector<double>&);
  bool tuple_predicate(const std::tuple<int,int,double>&,const std::tuple<int,int,double>&);
  bool pair_predicate_dbl(const std::pair<int,double>&,const std::pair<int,double>&);
  bool pair_predicate_int(const std::pair<int,int>&,const std::pair<int,int>&);
  int element(const std::set<int>&);
  bool next_combination(std::vector<int>&,int);
  void complement(const std::set<int>&,std::set<int>&,int,int);
  int coincidence(const std::set<int>&,const std::set<int>&);
  bool double_equality(double,double);
  void trim(std::string&);
  void set_wcomponent_values(int,bool);
  bool bfs(const pair_index&,int,int,int,int*);
  int network_flow(pair_index&,int,int,int);

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
    if (x == y) throw std::invalid_argument("The key arguments must be distinct integers!");

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
    if (v.empty()) throw std::invalid_argument("The key argument must not be an empty vector!");

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
    if (S.empty()) throw std::invalid_argument("The key argument must not be an empty set!");

    unsigned int n = S.size();
    std::stringstream s;
    std::set<int>::const_iterator it;

    if (n == 1) {
      it = S.begin();
      s << *it;
    }
    else {
      unsigned int i = 0;
      for(it=S.begin(); it!=S.end(); ++it) {
        if (i < (n-1)) {
          s << *it << ":";
        }
        else {
          s << *it;
        }
        ++i;
      }
    }
    return s.str();
  }
}
#endif
