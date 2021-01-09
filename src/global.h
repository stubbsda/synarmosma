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
  /// Standard LAPACK function to solve a general double precision linear system \f$Ax = b\f$. 
  void dgemm_(char*,char*,int*,int*,int*,double*,double*,int*,double*,int*,double*,double*,int*);
  /// Standard LAPACK function to compute the eigenvalues and eigenvectors of a symmetric double precision matrix.
  void dsyev_(char*,char*,int*,double*,int*,double*,double*,int*,int*);
  /// Standard LAPACK function to compute the eigenvalues and eigenvectors of a general double precision matrix.
  void dgesv_(int*,int*,double*,int*,int*,double*,int*,int*);
  /// Standard LAPACK function to compute the eigenvalues and eigenvectors of a general double precision complex matrix.
  void zgesv_(int*,int*,std::complex<double>*,int*,int*,std::complex<double>*,int*,int*);
  /// Standard LAPACK function to compute the inverse of a general double precision matrix using LU factorization.  
  void dgetri_(int*,double*,int*,int*,double*,int*,int*);
}

/** @file */

namespace SYNARMOSMA {
  /// A Boost type which is used to provide an indexed table for \f$n\f$-simplices for \f$n\ge 1\f$.
  using hash_map = boost::unordered_map<std::set<int>,int>;
  /// A Boost type used to store an indexed table of edges (each of which is std::pair<int,int>) for a graph. 
  using pair_index = boost::unordered_map<std::pair<int,int>,int>;

#ifdef __LP64__
  /// A portable signed 64 bit integer.
  using INT64 = signed long;
  /// A portable unsigned 64 bit integer.
  using UINT64 = unsigned long;
#else
  /// A portable signed 64 bit integer.
  using INT64 = signed long long;
  /// A portable unsigned 64 bit integer.
  using UINT64 = unsigned long long;
#endif

#ifdef DISCRETE
  /// This constant represents the smallest possible energy value for a vertex.
  const double energy_quantum = 1e-7;
  /// This constant represents the smallest possible spatial separation between two vertices.
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
  /// This enumerated class lists the three relational states that may exist among two objects x and y, namely 
  /// x lies before y, x lies after y or x and y are in some manner incommensurate and thus no relation may be 
  /// said to connect them.  
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
  /// This function returns the binomial coefficient, i.e. the coefficient of the \f$x^k\f$ term in the expansion of \f$(1+x)^n\f$, where \f$n\f$ and \f$k\f$ are the first and second arguments respectively.
  UINT64 binomial(int,int);
  /// This function computes the next combination of the first n integers taken k at a time, where k is the length of its first argument and n is the second argument. This next combination is stored in the first argument when the return value is true; if the return value is false it means that the combinations have been exhausted.
  bool next_combination(std::vector<int>&,int);
  /// This function generates all the combinations of the elements of its first argument, taken r (the second argument) at a time, storing them all in the final argument.
  int combinations(const std::set<int>&,int,std::vector<int>&);
  /// This recursive function computes all \f$n!\f$ permutations of the \f$n\f$ elements of its second argument, appending each such permutation to the first argument. The third argument is a remainder vector and the fourth the value of \f$n\f$.
    // permutation to the array "output".
  void permute(std::vector<std::vector<int> >&,std::vector<int>&,const std::vector<int>&,int);
  /// This function computes the induced orientation on the \f$n+1\f$ faces of a \f$n\f$-simplex. The first argument is \f$n+1\f$, the second argument is the set of vertices of the \f$n\f$-simplex, the third argument the parity of this same simplex, the fourth argument an index table that maps each face to an element of the final argument, the vector storing the induced orientation on each face. 
  void induced_orientation(int,const std::vector<int>&,int,const hash_map&,int*);
  /// This function takes an array of vectors of length equal to the first argument and computes which of these vectors is a neighbour of the vector specified by the third argument (an index of the array). The index of these neighbours is then stored in the final argument of the function. 
  void assemble_neighbours(int,const std::vector<int>*,int,std::vector<int>&);
  /// This function accepts two integer vectors of equal length, counts the number N of elements which are not the same and then returns +1 if N/2 is even and -1 if N/2 is odd. 
  int parity(const std::vector<int>&,const std::vector<int>&);
  /// This function assumes its argument is an integer vector of length \f$N\f$ consisting of a permutation of the set \f$\{0,1,\dots,N-1\}\f$ and then returns the parity of this permutation.
  int parity(const std::vector<int>&);
  /// This function establishes an ordering on the pairs \f$(n,x)\f$ and \f$(m,y)\f$ using the criterion \f$ x < y\f$, i.e. the second (floating point) element. 
  bool pair_predicate_dbl(const std::pair<int,double>&,const std::pair<int,double>&);
  /// This function establishes an ordering on the integral pairs \f$(n,m)\f$ and \f$(j,k)\f$ using the criterion \f$ m < k\f$, i.e. the second element. 
  bool pair_predicate_int(const std::pair<int,int>&,const std::pair<int,int>&);
  /// This function establishes an ordering on the triples \f$(n,m,x)\f$ and \f$(j,k,y)\f$ using the criterion \f$ x < y\f$, i.e. the third (floating point) element. 
  bool tuple_predicate(const std::tuple<int,int,double>&,const std::tuple<int,int,double>&);
  /// This function computes the prime factorization of its first argument, storing the result of as a vector of pairs representing the prime number and its power.
  void factorize(long,std::vector<std::pair<long,int> >&);
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
  /// This function checks that its argument is a singleton set and if so, returns the set's unique element.
  int element(const std::set<int>&);
  /// This function loops over the integers from 0 to N - 1 where N is the final argument, skipping the integer equal to the function's third argument. For each such integer, if the function doesn't find it in its first argument it is inserted in the second argument. 
  void complement(const std::set<int>&,std::set<int>&,int,int);
  /// This function returns the incidence number of its two arguments, assumed to be of length N and N+1 respectively. If the first argument is obtained from the second by dropping the \f$k\f$-th element, the function returns \f$(-1)^k\f$. If there is at least one element of the first argument not in the second, the function returns zero.
  int coincidence(const std::set<int>&,const std::set<int>&);
  /// This function determines if its two arguments differ by less than machine epsilon, returning true when this is the case and false otherwise.
  bool double_equality(double,double);
  /// This function uses a breadth first search to determine if a given source vertex (second argument) can reach a given sink vertex (third argument) in a graph with a table of edges (first argument) and a given number of vertices (fourth argument). If so the function returns true and the final argument lists the path of vertices from the source to the sink. 
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
