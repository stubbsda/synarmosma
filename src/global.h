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
#include <chrono>
#include <complex>
#include <algorithm>
#include <ctime>
#include <string>
#include <limits>
#include <tuple>
#include <queue>
#include <sys/time.h>
// OpenMP header file if applicable...
#ifdef _OPENMP
#include <omp.h>
#endif
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

/** @file */

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

  /// This function returns \f$b^n\f$, where \f$b > 0\f$ and \f$n > 0\f$ are its first and second arguments, respectively.
  UINT64 ipow(int,int);
  /// This function returns the factorial of its argument, which must be non-negative.
  UINT64 factorial(int);
  UINT64 binomial(int,int);
  int combinations(const std::set<int>&,int,std::vector<int>&);
  void permute(std::vector<std::vector<int> >&,std::vector<int>&,const std::vector<int>&,int);
  void assemble_neighbours(int,const std::vector<int>*,int,std::vector<int>&);
  void induced_orientation(int,const std::vector<int>&,int,const hash_map&,int*);
  void factorize(long,std::vector<std::pair<long,int> >&);
  int parity(const std::vector<int>&,const std::vector<int>&);
  int parity(const std::vector<int>&);
  void vertex_difference(int,int,std::vector<double>&);
  /// This function returns the arithmetic mean of the elements of its argument. 
  double arithmetic_mean(const std::vector<double>&);
  /// This function removes all of the leading and trailing white space from its argument.
  void trim(std::string&);
  /// This function tokenizes its first argument using the second argument as the delimiter, storing the tokens in the final argument.
  void split(const std::string&,char,std::vector<std::string>&);
  /// This function takes its vector argument of length n and creates n sequences of length (n-1) obtained by successively removing one element from the vector, stored in the second argument.
  void split(const std::vector<int>&,std::vector<int>&);
  /// This function adds to or subtracts from its first argument \f$x\f$ the quantity \f$L\f$ as many times as needed so that \f$x\in [0,L]\f$.
  double renormalize(double,double);
  /// This function writes an instance of the NTL::ZZ class (the second argument) to a binary output stream (first argument), returning the number of bytes written.
  int write_ZZ(std::ofstream&,const NTL::ZZ&);
  /// This function reads an instance of the NTL::ZZ class (the second argument) from a binary input stream (first argument), returning the number of bytes read.
  int read_ZZ(std::ifstream&,NTL::ZZ&);
  /// This function returns the square root of the sum of the squares of the elements of its first argument, an array whose length is equal to the second argument.
  double norm(const double*,int);
  /// This function returns the square root of the sum of the squares of the elements of its argument.
  double norm(const std::vector<double>&);
  /// This function computes the cross product of its two first arguments, a pair of 3-vectors, and writes the answer to the function's final argument.
  void cross_product(const std::vector<double>&,const std::vector<double>&,std::vector<double>&);
  bool tuple_predicate(const std::tuple<int,int,double>&,const std::tuple<int,int,double>&);
  bool pair_predicate_dbl(const std::pair<int,double>&,const std::pair<int,double>&);
  bool pair_predicate_int(const std::pair<int,int>&,const std::pair<int,int>&);
  /// This function checks that its argument is a singleton set and if so, returns the set's unique element.
  int element(const std::set<int>&);
  bool next_combination(std::vector<int>&,int);
  void complement(const std::set<int>&,std::set<int>&,int,int);
  int coincidence(const std::set<int>&,const std::set<int>&);
  /// This function determines if its two arguments differ by less than machine epsilon, returning true when this is the case and false otherwise.
  bool double_equality(double,double);
  bool breadth_first_search(const pair_index&,int,int,int,int*);
  /// This function uses the Ford-Fulkerson algorithm to return the maximum flow between two vertices of a graph. The function makes use of the breadth_first_search() function while the two vertices are the second and third arguments. The first argument is the index of edge capacities and the final argument the number of vertices in the graph. 
  int network_flow(const pair_index&,int,int,int);
  /// This function converts its first argument \f$0\le \rho\le \pi/2\f$ into three RGB colour intensities, defined as \f$255\sin^2\rho\f$, \f$255\sin2\rho\f$ and \f$255\cos^2\rho\f$ respectively.  
  inline void RGB_intensity(double rho,unsigned char* output)  
  {
    if (rho < -std::numeric_limits<double>::epsilon() || (rho - M_PI/2.0) > std::numeric_limits<double>::epsilon()) throw std::invalid_argument("Illegal angular argument in SYNARMOSMA::RGB_intensity!");
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
  /// This function returns its argument, a complex number \f$ z = x + iy\f$, as a string formatted "(x,y)". 
  inline std::string complex2string(const std::complex<double>& z)
  {
    std::string output = "(" + std::to_string(z.real()) + "," + std::to_string(z.imag()) + ")";
    return output; 
  }
  /// This function returns its integer argument as a string.
  inline std::string make_key(int x)
  {
    std::stringstream s;
    s << x;
    return s.str();
  }
  /// This function returns its two distinct arguments x and y as a colon-separated string, "x:y".
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
  /// This function returns a string which is a colon-separated list of the elements of the vector that is its argument.
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
  /// This function returns a string which is a colon-separated list of the elements of the set that is its argument.
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
