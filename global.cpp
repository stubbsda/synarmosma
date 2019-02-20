#include "random.h"

SYNARMOSMA::Random RND;

namespace SYNARMOSMA {
  bool tuple_predicate(const std::tuple<int,int,double>& lhs,const std::tuple<int,int,double>& rhs)
  {
    return std::get<2>(lhs) < std::get<2>(rhs);
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
    if (r > n) throw std::invalid_argument("The combination size must not exceed the number of elements!");

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

  int write_ZZ(std::ofstream& s,const NTL::ZZ& p)
  {
    int i,count = 0,n = NTL::NumBytes(p);
    bool neg = (p < 0) ? true : false;
    unsigned char q[n];

    s.write((char*)(&neg),sizeof(bool)); count += sizeof(bool);
    s.write((char*)(&n),sizeof(int)); count += sizeof(int);   
    NTL::BytesFromZZ(q,p,n);
    for(i=0; i<n; ++i) {
      s.write((char*)(&q[i]),sizeof(char)); count += sizeof(char);
    } 
    return count;
  }

  int read_ZZ(std::ifstream& s,NTL::ZZ& p)
  {
    int i,n,count = 0;
    bool neg;

    s.read((char*)(&neg),sizeof(bool)); count += sizeof(bool);
    s.read((char*)(&n),sizeof(int)); count += sizeof(int);
    unsigned char q[n];
    for(i=0; i<n; ++i) {
      s.read((char*)(&q[i]),sizeof(char)); count += sizeof(char);
    }
    p = NTL::ZZFromBytes(q,n);
    if (neg) p = -p;

    return count;
  }

  int parity(const std::vector<int>& v1,const std::vector<int>& v2)
  {
    if (v1.size() != v2.size()) throw std::invalid_argument("The parity arguments must be vectors of the same length!");

    unsigned int i,n = 0;
    for(i=0; i<v1.size(); ++i) {
      if (v1[i] != v2[i]) n++;
    }
    if ((n % 2) != 0) throw std::runtime_error("The raw parity distance is not an even number!");
    n /= 2;
    int output = (n%2 == 0) ? 1 : -1;
    return output;
  }

  void cross_product(const std::vector<double>& x,const std::vector<double>& y,std::vector<double>& z)
  {
    if (x.size() != 3 || y.size() != 3) throw std::invalid_argument("The cross-product is only defined for 3-vectors!");

    z.clear();
    z.push_back(0.0); z.push_back(0.0); z.push_back(0.0);
    z[0] = x[1]*y[2] - x[2]*y[1];
    z[1] = x[2]*y[0] - x[0]*y[2]; 
    z[2] = x[0]*y[1] - x[1]*y[0];
  }

  void assemble_neighbours(int nfacets,const std::vector<int>* facets,int in1,std::vector<int>& N) 
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
      p = ((i+1)%2 == 0) ? 1 : -1; 
      // n is the index of vx in the array of faces...
      qt = face_index.find(vx);
      if (qt == face_index.end()) throw std::runtime_error("The simplex face does not exist in the index table!");
      n = qt->second;
      output[n] = parity*p;
      vx.clear();
    }
  }
  /*
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
  */
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

  int parity(const std::vector<int>& symbols)
  {
    // A function to compute the parity of a permutation of (0,1,...,n-1) by tracing 
    // the trajectory of each element in the permutation and computing its length  
    const int n = (signed) symbols.size();
    int i,cpoint = 0,p = 1,length = 0,visited[n];

    for(i=0; i<n; ++i) {
      visited[i] = 0;
    }
    do {
      for(i=0; i<n; ++i) {
        if (symbols[i] == cpoint) {
          if (visited[i] == 0) {
            visited[i] = 1;
            cpoint = i;
            length++;
          }
          else {
            cpoint = -1;
          }
          break;
        }
      }
      if (cpoint == -1) {
        p *= ((length-1)%2 == 0) ? 1 : -1;
        length = 0;
        for(i=0; i<n; ++i) {
          if (visited[i] == 0) cpoint = i;
        }
        if (cpoint == -1) break;
      }
    } while(true);

    return p;
  }

  void permute(std::vector<std::vector<int> >& output,std::vector<int>& symbols,const std::vector<int>& rem,int n)
  {
    // A method that computes all n! permutations of the elements of the n-vector "symbols", appending each 
    // permutation to the array "output".
    if (rem.size() == 0) {
      output.push_back(symbols);
    }
    else {
      int i,j;
      std::vector<int> temp;
      const int m = (signed) rem.size();

      for(i=0; i<m; ++i) {
        symbols.push_back(rem[i]);
        for(j=0; j<m; ++j) {
          if (j == i) continue;
          temp.push_back(rem[j]);
        }
        permute(output,symbols,temp,n);
        temp.clear();
        symbols.pop_back();
      }
    }
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

  UINT64 binomial(int n,int k)
  {
    int i,ulimit = (k < (n-k)) ? k : n - k;
    UINT64 output,num = 1,den = 1;

    for(i=1; i<=ulimit; ++i) {
      num *= n + 1 - i;
      den *= i;
    }
    output = num/den;
    return output;
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
    if (v2.size() != (1 + v1.size())) throw std::invalid_argument("The coincidence arguments must differ in length by one!");
    int i,in1,output;
    bool found = true;
    std::set<int>::const_iterator it,jt;

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

  bool bfs(const pair_index& rgraph,int source,int sink,int nvertex,int* parent)
  {
    int u,v,kappa;
    std::pair<int,int> pr;
    pair_index::const_iterator pt;

    // Create a visited array and mark all vertices as not visited
    bool visited[nvertex];
    for(v=0; v<nvertex; ++v) {
      visited[v] = false;
    } 
    // Create a queue, enqueue source vertex and mark source vertex
    // as visited
    std::queue<int> q;
    q.push(source);
    visited[source] = true;
    parent[source] = -1; 
    // Standard BFS Loop
    do {
      u = q.front();
      q.pop();
      for(v=0; v<nvertex; ++v) {
        pr.first = u;
        pr.second = v;
        pt = rgraph.find(pr);
        kappa = (pt == rgraph.end()) ? 0 : pt->second;
        if (!visited[v] && kappa > 0) {
          q.push(v);
          parent[v] = u;
          visited[v] = true;
        }
      }
    } while(!q.empty()); 
    // If we reached sink in BFS starting from source, then return
    // true, else false
    return(visited[sink] == true);
  }

  int network_flow(pair_index& rgraph,int source,int sink,int nvertex)
  {
    // An implementation of the Ford-Fulkerson algorithm with breadth-first search 
    // (the Edmonds-Karp variation) for computing the maximum flow in a network with 
    // integer capacity and flow
    int u,v,g,path_flow,parent[nvertex],max_flow = 0;
    std::pair<int,int> pr;
    pair_index::const_iterator pt;

    // Augment the flow while there is path from source to sink
    while(bfs(rgraph,source,sink,nvertex,parent)) {
      // Find minimum residual capacity of the edges along the
      // path filled by BFS or we can say find the maximum flow
      // through the path found.
      path_flow = std::numeric_limits<int>::max();
      for(v=sink; v!=source; v=parent[v]) {
        u = parent[v];
        pr.first = u; pr.second = v;
        pt = rgraph.find(pr);
        g = (pt == rgraph.end()) ? 0 : pt->second;
        path_flow = std::min(path_flow,g);
      }
      // Update residual capacities of the edges and reverse edges
      // along the path
      for(v=sink; v!=source; v=parent[v]) {
        u = parent[v];
        pr.first = u; pr.second = v;
        rgraph[pr] -= path_flow;
        pr.first = v; pr.second = u;
        rgraph[pr] += path_flow;
      }
      // Add path flow to overall flow
      max_flow += path_flow;
    } 
    // Return the overall flow
    return max_flow;
  }
}
