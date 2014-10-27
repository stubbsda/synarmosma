#include "global.h"

SYNARMOSMA::Random RND;

namespace SYNARMOSMA {
  bool tuple_predicate(const boost::tuple<int,int,double>& lhs,const boost::tuple<int,int,double>& rhs)
  {
    return lhs.get<2>() < rhs.get<2>();
  }

  bool pair_predicate_dbl(const std::pair<int,double>& lhs,const std::pair<int,double>& rhs)
  {
    return lhs.second < rhs.second;
  }

  bool pair_predicate_int(const std::pair<int,int>& lhs,const std::pair<int,int>& rhs)
  {
    return lhs.second < rhs.second;
  }

  bool double_equality(double x,double y)
  {
    return std::abs(x - y) <= ( (std::abs(x) < std::abs(y) ? std::abs(y) : std::abs(x))*std::numeric_limits<double>::epsilon());
  }

  int combinations(const std::set<int>& S,int r,std::vector<int>& C)
  {
    // A function to generate all the combinations of the elements of S, taken 
    // r at a time.
    const int n = (signed) S.size();
    assert(r <= n);

    int i;
    std::vector<int> v;
    std::set<int>::const_iterator it;

    C.clear();

    for(i=0; i<r; ++i) {
      v.push_back(1);
    }
    for(i=r; i<n; ++i) {
      v.push_back(0);
    }
    do {
      i = 0;
      for (it=S.begin(); it!=S.end(); ++it) {
        if (v[i] == 1) C.push_back(*it);
        ++i;
      }
    } while (std::next_permutation(v.begin(),v.end()));
    i = (signed) C.size()/r;
    return i;
  }

  void factorize(long n,std::vector<std::pair<long,int> >& factors)
  {
    int i,d,np;
    bool done;
    long q,f,p = 1,L = n/2;
    std::vector<long> psequence;

    factors.clear();

    do {
      q = NTL::NextPrime(p);
      psequence.push_back(q);
      p = q+1;
      if (p >= L) break;
    } while(true);
    np = (signed) psequence.size();

    for(i=0; i<np; ++i) {
      f = psequence[i];
      d = 0;
      done = false;
      do {
        if ((n % f) == 0) { 
          n = n/f;
          d++;
        }
        else {
          done = true;
        }
      } while(!done);
      if (d > 0) factors.push_back(std::pair<long,int>(f,d));
      if (n == 1) break;
   }
  } 

  int parity(const std::vector<int>& v1,const std::vector<int>& v2)
  {
    assert(v1.size() == v2.size());
    int i,n = 0;
    for(i=0; i<(signed) v1.size(); ++i) {
      if (v1[i] != v2[i]) n++;
    }
    assert((n % 2) == 0);
    n /= 2;
    return ipow(-1,n);
  }

  void cross_product(const std::vector<double>& x,const std::vector<double>& y,std::vector<double>& z)
  {
    assert(x.size() == 3);
    assert(y.size() == 3);
    z.clear();
    z.push_back(0.0); z.push_back(0.0); z.push_back(0.0);
    z[0] = x[1]*y[2] - x[2]*y[1];
    z[1] = x[2]*y[0] - x[0]*y[2]; 
    z[2] = x[0]*y[1] - x[1]*y[0];
  }

  void get_neighbours(int nfacets,const std::vector<int>* facets,int in1,std::vector<int>& N) 
  {
    std::vector<int> base = facets[in1];
    int i,j,k,nfound,D = (signed) base.size();
    bool found;

    N.clear();
    for(i=0; i<nfacets; ++i) {
      if (i == in1) continue;
      // If facets[i] and base differ by a single vertex, they are neighbours 
      nfound = 0;
      for(j=0; j<D; ++j) {
        found = false;
        for(k=0; k<D; ++k) {
          if (base[j] == facets[i][k]) {
            found = true;
            break;
          }
        }
        if (found) nfound++;
      }
      if (nfound < (D-1)) continue;
      N.push_back(i);
    }
  }

  void induced_orientation(int nf,const std::vector<int>& S,int parity,const hash_map& face_index,int* output) 
  {
    // Break out the D+1 faces of the D-simplex S along with their parity, 
    int i,j,p,n,D = (signed) S.size();
    hash_map::const_iterator qt;
    std::set<int> vx;

    for(i=0; i<nf; ++i) {
      output[i] = 0;
    }
    for(i=0; i<D; ++i) {
      for(j=0; j<D; ++j) {
        if (j == i) continue;
        vx.insert(S[j]);
      }
      p = ipow(-1,i+1);
      // n is the index of vx in the array of faces...
      qt = face_index.find(vx);
      assert(qt != face_index.end());
      n = qt->second;
      output[n] = parity*p;
      vx.clear();
    }
  }

  void convert(unsigned char* p,int n)
  {
    unsigned int i;
    unsigned char* temp = (unsigned char*)(&n);
    for(i=0; i<sizeof(int); ++i) {
      p[i] = temp[i];
    }
  }

  void convert(unsigned char* p,float x)
  {
    unsigned int i;
    unsigned char* temp = (unsigned char*)(&x);
    for(i=0; i<sizeof(float); ++i) {
      p[i] = temp[i];
    }
  }

  void split(const std::string& s,char delim,std::vector<std::string>& elements) 
  {
    std::stringstream ss(s);
    std::string item;

    elements.clear();

    while(std::getline(ss,item,delim)) {
      elements.push_back(item);
    }
  }

  int element(const std::set<int>& vx)
  {
    if (vx.empty()) return -1;
    if (vx.size() > 1) return -1;
    int out = -1;
    std::set<int>::const_iterator it;
    for(it=vx.begin(); it!=vx.end(); ++it) {
      out = *it;
    }
    return out;
  }

  bool next_combination(std::vector<int>& vx,int n)
  {
    int i,j;
    const int k = (signed) vx.size();
    vx[k-1] += 1;
    for(i=k-1; i>0; --i) {
      if (vx[i] >= (n-k+1+i)) {
        vx[i-1] += 1;
      }
      else {
        break;
      }
    }

    if (vx[0] > (n - k)) return false;

    for(j=i+1; j<k; ++j) {
      vx[j] = vx[j-1] + 1;
    }
    return true;
  }

  double norm(const double* x,int n)
  {
    double out = 0.0;
    for(int i=0; i<n; ++i) {
      out += x[i]*x[i];
    }
    return std::sqrt(out);
  }

  double norm(const std::vector<double>& x)
  {
    double alpha,output = 0.0;
    std::vector<double>::const_iterator vit;
    for(vit=x.begin(); vit!=x.end(); ++vit) {
      alpha = *vit;
      output += alpha*alpha;
    }
    return std::sqrt(output);
  }

  double binomial(int n,int k)
  {
    int i,f,output = n;  
    for(i=1; i<k; ++i) {
      output *= (n - i);
    }
    f = 1;
    for(i=2; i<=k; ++i) {
      f *= i;
    }
    output /= f;
    return double(output);
  }

  UINT64 factorial(int x)
  {
    UINT64 output = 1;
    for(int i=2; i<=x; ++i) {
      output *= (UINT64) i;
    }
    return output;
  }

  double dmap(double x,double L)
  {
    // A function to map the variable x to the range [0,L], by 
    // adding or subtracting multiples of L
    double output = x;
    if (output < 0.0) {
      do {
        output += L;
        if (output > 0.0) break;
      } while(true);
    }
    else if (output > L) {
      do {
        output -= L;
        if (output < L) break;
      } while(true);
    }
    return output;
  }

  double arithmetic_mean(const std::vector<double>& input)
  {
    int i,n = (signed) input.size();
    double sum = 0.0;
    for(i=0; i<n; ++i) {
      sum += input[i];
    }
    return sum/double(n);
  }

  UINT64 ipow(int b,int n)
  {
    UINT64 output = (UINT64) b;
    if (n == 0) return 1;
    for(int i=1; i<n; ++i) {
      output *= (UINT64) b;
    }
    return output;
  }

  void trim(std::string& str)
  {
    std::string::size_type pos = str.find_last_not_of(' ');
    if (pos != std::string::npos) {
      str.erase(pos+1);
      pos = str.find_first_not_of(' ');
      if (pos != std::string::npos) str.erase(0,pos);
    }
    else {
      str.erase(str.begin(),str.end());
    }
  }

  void split(const std::vector<int>& p,std::vector<int>& c)
  {
    // This function takes the integer vector p of length n and
    // creates n sequences of length (n-1) obtained by removing
    // one element from p.
    int i,j,n = (signed) p.size();
    for(i=0; i<n; ++i) {
      for(j=0; j<n; ++j) {
        if (i == j) continue;
        c.push_back(p[j]);
      }
    }
  }

  void invert(const double* A,double* B,int n)
  {
    // Invert the size*size matrix A into B
    int i,info,lwork = n*n;
    int ipiv[n+1];
    double work[lwork];

    for(i=0; i<n*n; ++i) {
      B[i] = A[i];
    }
    dgetrf_(&n,&n,B,&n,ipiv,&info);
    dgetri_(&n,B,&n,ipiv,work,&lwork,&info);
    assert(info == 0); 
  }

  double determinant(const double* A,int n) 
  {
    // A function that calculates the determinant of the n x n matrix A, 
    // using LAPACK to first compute an LU factorization of A.
    int i,info,pivots[n];
    bool permuted = false;
    double B[n*n],output = 1.0;

    for(i=0; i<n*n; ++i) {
      B[i] = A[i];
    }
    dgetrf_(&n,&n,B,&n,pivots,&info);
    assert(info >= 0);
    for(i=0; i<n; ++i) {
      output *= B[n*i+i];
      if (pivots[i] != (1+i)) permuted = !permuted;
    }
    output = (permuted) ? -output : output;
    return output;  
  }

  int coincidence(const std::set<int>& v1,const std::set<int>& v2)
  {
    // The incidence number concerns the relation between two
    // lists of integers, v2 of length n and v1 of length n-1.
    // If there is at least one element of v1 not found in v2,
    // the incidence number is zero. If we define the longer
    // list by
    // v2 = (v_1,v_2,...,v_(1+n))
    // and the shorter by
    // v1 = (v_1,v_2,...,v_k,v_(k+1),...v_(1+n))
    // that is, v1 is obtained from v2 by dropping the k-th
    // element of v2, then the incidence number is (-1)^k.
    int i,in1,output;
    std::set<int>::const_iterator it,jt;
    bool found = true;

    assert((1+v1.size()) == v2.size());

    for(it=v1.begin(); it!=v1.end(); ++it) {
      in1 = *it;
      jt = std::find(v2.begin(),v2.end(),in1);
      if (jt == v2.end()) {
        found = false;
        break;
      }
    }
    if (!found) return 0;

    i = 1;
    for(it=v2.begin(); it!=v2.end(); ++it) {
      in1 = *it;
      jt = std::find(v1.begin(),v1.end(),in1);
      if (jt == v1.end()) {
        output = (i%2) ? 1 : -1;
        return output;
      }
      ++i;
    }
    // We should never get here in fact...
    return 0;
  }

  void complement(const std::set<int>& input,std::set<int>& output,int self,int modulus)
  {
    std::set<int>::const_iterator it;
    output.clear();
    for(int i=0; i<self; ++i) {
      it = std::find(input.begin(),input.end(),i);
      if (it == input.end()) output.insert(i);
    }
    for(int i=self+1; i<modulus; ++i) {
      it = std::find(input.begin(),input.end(),i);
      if (it == input.end()) output.insert(i);
    }
  }
}
