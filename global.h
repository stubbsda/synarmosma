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
  using base_generator_type = boost::mt19937;
  using hash_map = boost::unordered_map<std::set<int>,int>;
  using edge_hash = boost::unordered_map<std::pair<int,int>,int>;

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
  class Monomial {
   public:
    kind coefficient;
    std::vector<std::pair<unsigned int,unsigned int> > exponents; 
  };

  enum class Relation
  {
      before,
      after,
      disparate
  };

  UINT64 ipow(int,int);
  UINT64 factorial(int);
  UINT64 binomial(int,int);
  int combinations(const std::set<int>&,int,std::vector<int>&);
  void permute(std::vector<std::vector<int> >&,std::vector<int>&,const std::vector<int>&,int);
  void get_neighbours(int,const std::vector<int>*,int,std::vector<int>&);
  void induced_orientation(int,const std::vector<int>&,int,const hash_map&,int*);
  void factorize(long,std::vector<std::pair<long,int> >&);
  int parity(const std::vector<int>&,const std::vector<int>&);
  int parity(const std::vector<int>&);
  void convert(unsigned char*,int);
  void convert(unsigned char*,float);
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
  bool bfs(const edge_hash&,int,int,int,int*);
  int network_flow(edge_hash&,int,int,int);

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
#ifdef DEBUG
    assert(x != y);
#endif
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
#ifdef DEBUG
    assert(!v.empty());
#endif
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
#ifdef DEBUG
    assert(!S.empty());
#endif
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
    bool poisson_allocated;
   public:
    Random();
    ~Random();
    inline void set_seed(int x) {assert(x >= 0); s = x; BGT.seed(s);};
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
    void generate_random_vector(std::vector<double>&,int,double,double);
    void generate_random_vector(std::vector<std::complex<double> >&,int,double,double,bool = false);
  };
}
#endif
