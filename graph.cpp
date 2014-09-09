#include "graph.h"

extern Random RND;

Graph::Graph() : Schema()
{
  nvertex = 0;
  nedge = 0;
}

Graph::Graph(int n) : Schema(n)
{
  // The complete graph on n vertices...
  int i,j;

  nedge = 0;
  for(i=0; i<n; ++i) {
    for(j=1+i; j<n; ++j) {
      add_edge(i,j);
    }
  }
}

Graph::Graph(int n,int c) : Schema(n)
{
  // Construct a scale-free graph with n vertices where 
  // the minimum degree is c
  int i,j,k,sum;
  double rho;

  assert(c >= 1);
 
  nedge = 0;
  for(i=0; i<c+1; ++i) {
    for(j=(signed) neighbours[i].size(); j<c; ++j) {
      do {
        k = RND.irandom(c+1);
        if (k == i) continue;
        if (neighbours[i].count(k) == 0) break;
      } while(true);
      add_edge(i,k);
    }
  }

  for(i=c+1; i<n; ++i) {
    for(j=0; j<c; ++j) {
      rho = RND.drandom();
      sum = 0;
      for(k=0; k<i; ++k) {
        sum += neighbours[k].size();
        if (rho < double(sum)/double(2*nedge)) break;
      }
      add_edge(i,k);
    }
  }
}

Graph::Graph(int n,double p) : Schema(n)
{
  // We will use the Erdős–Rényi random graph model (the G(n,p) variant) to
  // assemble a random graph, with n = initial_size...
  int i,j;
  double alpha;

  nedge = 0;
  for(i=0; i<n; ++i) {
    // Compute a lower bound for the valence of this vertex
    for(j=1+i; j<n; ++j) {
      alpha = RND.drandom();
      if (alpha > p) continue;
      add_edge(i,j);
    }
  }
}

Graph::Graph(const Graph& source)
{
  clear();
  nvertex = source.nvertex;
  neighbours = source.neighbours;
  nedge = source.nedge;
}

Graph::~Graph()
{

}

void Graph::core(Graph* G,int k) const
{
  // A method to compute the k-core of the current graph
  // The 1-core is just the graph itself...
  assert(k > 1);

  if (k > max_degree()) {
    // The k-core in this case is null
    G->clear();
    return;
  }
  else if (k < min_degree()) {
    // The k-core is just the graph as a whole in this case...
    G->nvertex = nvertex; 
    G->nedge = nedge;
    G->neighbours = neighbours;
    return;
  }

  int i;
  bool found;

  G->nvertex = nvertex; 
  G->nedge = nedge;
  G->neighbours = neighbours;

  do {
    if (G->max_degree() < k) {
      // The k-core will necessarily be null in this case, so we can exit
      G->clear();
      break;
    }   
    found = false;
    for(i=0; i<G->nvertex; ++i) {
      if ((signed) G->neighbours[i].size() < k) {
        G->amputation(i);
        found = true;
        break;
      }
    }
    if (!found) break;
  } while(true);
}

bool Graph::planar() const
{
  if (nvertex <= 2) return true;
  if (nedge > (3*nvertex - 6)) return false;
  // Now the hard case where we need to do some work
  // to get the answer...
  bool output = true;

  return output;
}

bool Graph::add_edge(int v,int u)
{
  if (u == v) return false;
  if (neighbours[v].count(u) == 0) {
    neighbours[v].insert(u);
    neighbours[u].insert(v);
    nedge++;
    return true;
  }
  return false;
}

void Graph::clear()
{
  nedge = 0;
  nvertex = 0;
  neighbours.clear();
}

bool Graph::amputation(int v)
{
  assert(v >= 0 && v < nvertex);
  int i,j;
  std::set<int> S = neighbours[v];
  std::set<int>::const_iterator it;

  for(it=S.begin(); it!=S.end(); ++it) {
    neighbours[*it].erase(v);
  }
  nedge -= S.size();
  for(i=0; i<nvertex; ++i) {
    if (i == v) continue;
    S.clear();
    for(it=neighbours[i].begin(); it!=neighbours[i].end(); ++it) {
      j = *it;
      if (j > v) {
        S.insert(j-1);
      }
      else {
        S.insert(j);
      }
    }
    neighbours[i] = S;
  }
  nvertex--;
  neighbours.erase(neighbours.begin() + v);
  return true;
}

bool Graph::fusion(int v,int u)
{
  if (u == v) return false;
  if (u == -1) {
    if (neighbours[v].empty()) return false;
    u = RND.irandom(neighbours[v]);
  }
  int i;
  std::set<int> S = neighbours[u];
  std::set<int>::const_iterator it;
  amputation(u);
  if (u < v) v = v - 1;
  for(it=S.begin(); it!=S.end(); ++it) {
    i = *it;
    if (u < i) i = i - 1;
    if (i == v) continue;
    add_edge(v,i);
  }
  return true;
}

bool Graph::foliation_x(int v,int u)
{
  if (u == v) return false;
  if (neighbours[v].empty()) return false;
  if (u == -1) u = RND.irandom(neighbours[v]);
  neighbours[v].erase(u);
  neighbours[u].erase(v);
  nedge--;
  return true;
}

bool Graph::foliation_m(int v,int u)
{
  if (u == -1) {
    do {
      u = RND.irandom(nvertex);
      if (u != v) break;
    } while(true);
  }
  return add_edge(v,u);
}

int Graph::fission_x(int v)
{
  int u = nvertex;
  nvertex++;
  add_edge(v,u);
  return u;
}

int Graph::fission_m(int v)
{
  int u = nvertex;
  std::set<int> S = neighbours[v];
  std::set<int>::const_iterator it;
  nvertex++;
  add_edge(v,u);
  for(it=S.begin(); it!=S.end(); ++it) {
    add_edge(u,*it);
  }
  return u;
}

void Graph::katz_centrality(std::vector<double>& output) const
{
  // A method to compute the Katz centrality of the vertices in the 
  // graph - we first need to compute the adjacency matrix and its 
  // largest eigenvalue, then we use a recursive process in order to 
  // compute the vector of centrality values.
  int i,info,nv = nvertex,nwork = 3*nvertex - 1;
  unsigned int j;
  char jtype = 'N';
  char uplo = 'U';
  double alpha,sum;
  double* AD = new double[nvertex*nvertex];
  double* w = new double[nvertex];
  double* work = new double[nwork];
  std::vector<double> x,xnew;
  std::vector<unsigned int>::const_iterator it;
  Binary_Matrix* A = new Binary_Matrix;
  const double beta = 1.0;

  compute_adjacency_matrix(A);
  for(i=0; i<nvertex*nvertex; ++i) {
    AD[i] = 0.0;
  }

  for(i=0; i<nvertex; ++i) {
    for(it=A->elements[i].begin(); it!=A->elements[i].end(); ++it) {
      j = *it;
      AD[nvertex*i+j] = 1.0;
    }
  }

  // If info is zero after this LAPACK call, then w will contain
  // the eigenvalues of A (which are real since A is symmetric and
  // real) in ascending order, thus w[nv-1] will be the largest
  dsyev_(&jtype,&uplo,&nv,AD,&nv,w,work,&nwork,&info);
  assert(info == 0);  

  for(i=0; i<nvertex; ++i) {
    x.push_back(1.0);
    xnew.push_back(0.0);
  }

  // We want alpha to be just under the threshold for convergence
  alpha = 0.9*(1.0/w[nvertex-1]);

  delete AD;
  delete w;
  delete work;

  do {
    for(i=0; i<nvertex; ++i) {
      sum = 0.0;
      for(it=A->elements[i].begin(); it!=A->elements[i].end(); ++it) {
        sum += x[*it];
      }
      xnew[i] = alpha*sum + beta;
    }
    // Test for convergence...
    sum = 0.0;
    for(i=0; i<nvertex; ++i) {
      sum += (xnew[i] - x[i])*(xnew[i] - x[i]);
    }
    if (sum < 0.0001) break;
    x = xnew;
  } while(true);

  delete A;  
  output = xnew;
}

double Graph::clustering_coefficient(int v) const
{
  // This method calculates the percentage of distinct pairs (u,w) 
  // of neighbours of v which are also connected directly
  if (neighbours[v].size() < 2) return 0.0;
  int i,j,n,np,nf = 0;
  std::vector<int> vx;
  std::set<int>::const_iterator it;

  for(it=neighbours[v].begin(); it!=neighbours[v].end(); ++it) {
    vx.push_back(*it);
  }
  n = (signed) vx.size();
  for(i=0; i<n; ++i) {
    for(j=1+i; j<n; ++j) {
      if (neighbours[vx[i]].count(vx[j]) > 0 && neighbours[vx[j]].count(vx[i]) > 0) nf++;
    }
  }
  np = n*(n-1)/2;
  return double(nf)/double(np);
}

bool Graph::biconnected() const
{
  // If there exists at least one vertex whose removal (along
  // with its incident edges) renders the graph disconnected,
  // return false, otherwise return true
  int i,j,k,offset[nvertex];
  std::set<int>::const_iterator it;
  Graph G;

  for(i=0; i<nvertex; ++i) {
    offset[i] = -1;
  }

  for(i=0; i<nvertex; ++i) {
    // Delete this vertex and all of its incident edges...
    for(j=0; j<nvertex; ++j) {
      if (i == j) continue;
      offset[j] = G.add_vertex();
    }
    for(j=0; j<nvertex; ++j) {
      if (i == j) continue;
      for(it=neighbours[j].begin(); it!=neighbours[j].end(); ++it) {
        k = *it;
        if (k == i) continue;
        G.add_edge(offset[j],offset[k]);
      }
    }
    // Then test the connectivity...
    if (!G.connected()) return false;
    G.clear();
  }
  return true;
}

double Graph::return_probability(int base,int length) const
{
  int i,j,h = 0,next,current;
  bool home;
  double rho;
  const int ntrials = 50;

  for(i=0; i<ntrials; ++i) {
    current = base;
    home = false;
    for(j=0; j<length; ++j) {
      next = RND.irandom(neighbours[current]);
      current = next;
      if (current == base) {
        home = true;
        break;
      }
    }
    if (home) h++;
  }

  rho = double(h)/double(ntrials);
  return rho;
}

void Graph::random_walk(double* mean,double* sdeviation,int D) const
{
  int i,j,v;
  double mu,rho = 0.0,sigma = 0.0;
  std::set<int> vx;
  const int ntrials = 15;
  const int L = int(std::pow(nvertex,1.0/double(D)));
  const int nbase = int(0.05*nvertex);
  double value[nbase];

  for(i=0; i<nbase; ++i) {
    do {
      v = RND.irandom(nvertex);
      if (vx.count(v) > 0) continue;
      vx.insert(v);
      break;
    } while(true);
    mu = 0.0;
    for(j=0; j<ntrials; ++j) {
      mu += return_probability(v,L);
    }
    mu = mu/double(ntrials);
    value[i] = mu;
    rho += mu;
  }

  rho = rho/double(nbase);
  for(i=0; i<nbase; ++i) {
    sigma = (value[i] - rho)*(value[i] - rho);
  }
  sigma = std::sqrt(sigma/double(nbase));
  *mean = rho;
  *sdeviation = sigma;
}

void Graph::degree_distribution(bool logarithmic,std::vector<double>& histogram) const
{
  assert(connected());
  int i;
  const int max = max_degree();
  int counter[1+max];

  for(i=0; i<=max; ++i) {
    counter[i] = 0;
  }

  for(i=0; i<nvertex; ++i) {
    counter[neighbours[i].size()] += 1;
  }

  histogram.clear();
  histogram.push_back(double(counter[1]));
  if (logarithmic) {
    // Use logarithmic binning, so intervals of size {1,2,4,8,16...}
    int lbound = 1,ubound = 2,sum,its = 1;
    double alpha;
    do {
      lbound *= 2;
      ubound *= 2;
      if (lbound > max) break;
      if (ubound > (1+max)) ubound = 1 + max;
      sum = 0;
      for(i=lbound; i<ubound; ++i) {
        sum += counter[i];
      }
      alpha = double(sum)/double(ubound - lbound);
      histogram.push_back(alpha);
      its++;
    } while(true);
  }
  else {
    // Use a uniform bin width of unity...
    for(i=2; i<=max; ++i) {
      histogram.push_back(double(counter[i]));
    }
  }
}

double Graph::percolation(bool site) const 
{
  assert(connected());
  int i,n,nc;
  double output;
  std::vector<int> csize,components;
  const double NE = double(nedge);
  const double NV = double(nvertex);
  
  Graph wcopy(*this);

  if (site) {
    // Site percolation - we remove vertices and their associated edges...
    do {
      n = RND.irandom(wcopy.nvertex);
      if (!wcopy.amputation(n)) continue;
      nc = wcopy.component_analysis(components);
      if (nc == 1) continue;
      // We need to see if the giant component still exists...
      for(i=0; i<nc; ++i) {
        csize.push_back(0);
      }
      for(i=0; i<wcopy.nvertex; ++i) {
        csize[components[i]] += 1;
      }
      // Sort the component sizes in ascending order
      std::sort(csize.begin(),csize.end());
      // The giant component has vanished if the largest component is less than 
      // twice the size of the next largest component...
      if (double(csize[nc-1])/double(csize[nc-2]) < 2.0) break;
      csize.clear();
    } while(true);
    output = double(wcopy.nvertex)/NV;
  }
  else {
    // Bond percolation - we only remove edges..
    do {
      n = RND.irandom(wcopy.nvertex);
      // Eliminate a random edge connected to vertex "n"...
      if (!wcopy.foliation_x(n,-1)) continue;
      nc = wcopy.component_analysis(components);
      if (nc == 1) continue;
      // We need to see if the giant component still exists...
      for(i=0; i<nc; ++i) {
        csize.push_back(0);
      }
      for(i=0; i<wcopy.nvertex; ++i) {
        csize[components[i]] += 1;
      }
      // Sort the component sizes in ascending order
      std::sort(csize.begin(),csize.end());
      // The giant component has vanished if the largest component is less than 
      // twice the size of the next largest component...
      if (double(csize[nc-1])/double(csize[nc-2]) < 2.0) break;
      csize.clear();
    } while(true);
    output = double(wcopy.nedge)/NE;
  }
  return output;
}

double Graph::cosine_similarity(int u,int v) const
{
  // The value for this is the (u,v) element of the square of the 
  // adjacency matrix divided by the square root of the product of 
  // degree(u) and degree(v)
  if (u == v) return 1.0;
  if (neighbours[u].empty() || neighbours[v].empty()) return 0.0;

  int i,sum = 0;
  Binary_Matrix* A = new Binary_Matrix;
  const double denominator = std::sqrt(double(neighbours[u].size()*neighbours[v].size()));
  compute_adjacency_matrix(A);
  // Now we need to square this binary matrix A...
  for(i=0; i<nvertex; ++i) {
    if (A->get(u,i) && A->get(i,v)) sum += 1;
  }
  delete A;
  return double(sum)/denominator;
}

double Graph::inverse_girth() const
{
  assert(nvertex > 0);
  int i,j,p,q,alpha,output = 1 + nedge,done[nvertex],parent[nvertex],dist[nvertex];
  std::set<int> current,next;
  std::set<int>::const_iterator it,jt;

  for(i=0; i<nvertex; ++i) {
    current.insert(i);
    parent[i] = -1;
    dist[i] = 0;
    for(j=0; j<nvertex; ++j) {
      done[j] = 0;
    }
    do {
      for(it=current.begin(); it!=current.end(); ++it) {
        p = *it;
        done[p] = 1;
        for(jt=neighbours[p].begin(); jt!=neighbours[p].end(); ++jt) {
          q = *jt;
          if (q == parent[p]) continue;
          if (done[q] == 1) {
            alpha = 1 + dist[p] + dist[q];
            if (alpha < output) output = alpha;
          }
          else {
            parent[q] = p;
            dist[q] = 1 + dist[p];
            next.insert(q);
          }
        }
      }
      if (next.empty()) break;
      current = next;
      next.clear();
    } while(true);
    current.clear();
  }
  double ginv = (output == (1+nedge)) ? 0.0 : 1.0/double(output);
  return ginv;
}

int Graph::cyclomatic_number() const
{
  return (nedge - nvertex + 1);
}

int Graph::max_degree() const
{
  int i,n,output = 0;
  for(i=0; i<nvertex; ++i) {
    n = (signed) neighbours[i].size();
    if (n > output) output = n;
  }
  return output;
}

int Graph::min_degree() const
{
  int i,n,output = 1 + nedge;
  for(i=0; i<nvertex; ++i) {
    n = (signed) neighbours[i].size();
    if (n < output) output = n;
  }
  return output;
}

double Graph::average_degree() const
{
  int i;
  double sum = 0.0;
  for(i=0; i<nvertex; ++i) {
    sum += double(neighbours[i].size());
  }
  return sum/double(nvertex);
}

double Graph::completeness() const
{
  double output = 0.0;
  if (nvertex == 1) return output;

  output = 2.0*double(nedge)/double(nvertex*(nvertex-1));
  return output;
}

void Graph::compute_adjacency_matrix(Binary_Matrix* A) const
{
  int i,j;
  std::set<int>::const_iterator it;

  A->initialize(nvertex,nvertex);
  for(i=0; i<nvertex; ++i) {
    for(it=neighbours[i].begin(); it!=neighbours[i].end(); ++it) {
      j = *it;
      if (j < i) continue;
      A->elements[i].push_back(j); 
      A->elements[j].push_back(i); 
    }
  }
}

double Graph::entwinement() const
{
  // This method produces a real number between 0 and 1 that measures the
  // degree of "labyrinthicity" of the graph
  if (nvertex == 1) return 0.0;

  assert(connected());

  int i,info,nv = nvertex,nwork = 3*nvertex - 1;
  unsigned int j;
  char jtype = 'N';
  char uplo = 'U';
  double output = 0.0;
  double* AD = new double[nvertex*nvertex];
  double* w = new double[nvertex];
  double* work = new double[nwork];
  std::vector<unsigned int>::const_iterator it;

  // First built the adjacency matrix
  Binary_Matrix* A = new Binary_Matrix;
  compute_adjacency_matrix(A);

  for(i=0; i<nvertex*nvertex; ++i) {
    AD[i] = 0.0;
  }

  for(i=0; i<nvertex; ++i) {
    for(it=A->elements[i].begin(); it!=A->elements[i].end(); ++it) {
      j = *it;
      AD[nvertex*i+j] = 1.0;
    }
  }
  delete A;
  // If info is zero after this LAPACK call, then w will contain
  // the eigenvalues of A (which are real since A is symmetric and
  // real) in ascending order, thus w[nv-1] will be the largest
  dsyev_(&jtype,&uplo,&nv,AD,&nv,w,work,&nwork,&info);

  // We know that w[nv-1] <= max_degree(), so we divide by this
  // to normalize the output
  if (info == 0) {
    double ds = sqrt(double(max_degree()));
    double da = average_degree();
    double m1 = (ds < da) ? da : ds;
    output = (w[nvertex-1] - m1)/(ds*ds - m1);
  }

  delete[] AD;
  delete[] w;
  delete[] work;

  return output;
}

double Graph::cyclic_resistance() const
{
  assert(connected());
  int i,j,info,nv = nvertex;
  int nwork = 5*nvertex;
  unsigned int l;
  double sum;
  std::set<int>::const_iterator it;
  std::vector<unsigned int>::const_iterator vt;
  Binary_Matrix* A = new Binary_Matrix;
  double* L = new double[nvertex*nvertex];
  double* C = new double[nvertex*nvertex];
  double* W = new double[nvertex*nvertex];
  double* work = new double[nwork];
  int* pivots = new int[nvertex];

  compute_adjacency_matrix(A);

  for(i=0; i<nvertex; ++i) {
    L[nvertex*i+i] = double(neighbours[i].size());
    for(j=i+1; j<nvertex; ++j) {
      L[nvertex*i+j] = 0.0;
      L[nvertex*j+i] = 0.0;
    }
  }
  for(i=0; i<nvertex; ++i) {
    for(it=neighbours[i].begin(); it!=neighbours[i].end(); ++it) {
      j = *it;
      if (j < i) continue;
      L[nvertex*i+j] = -1.0;
      L[nvertex*j+i] = -1.0;
      continue;
    }
  }
  for(i=0; i<nvertex*nvertex; ++i) {
    C[i] = L[i];
  }
  dgetri_(&nv,C,&nv,pivots,work,&nwork,&info);

  for(i=0; i<nvertex; ++i) {
    W[nvertex*i+i] = 0.0;
    for(j=i+1; j<nvertex; ++j) {
      sum = 1.0/(C[nvertex*i+i] + C[nvertex*j+j] - (C[nvertex*i+j] + C[nvertex*j+i]));
      W[nvertex*i+j] = sum;
      W[nvertex*j+i] = sum;
    }
  }
  for(i=0; i<nvertex; ++i) {
    for(j=0; j<nvertex; ++j) {
      sum = 0.0;
      for(vt=A->elements[j].begin(); vt!=A->elements[j].end(); vt++) {
        l = *vt;
        sum += W[nvertex*l+i];
      }
      L[nvertex*i+j] = sum;
    }
  }
  sum = 0.0;
  for(i=0; i<nvertex; ++i) {
    sum += L[nvertex*i+i];
  }
  sum -= nedge;

  delete[] work;
  delete[] pivots;
  delete[] L;
  delete[] C;
  delete[] W;
  delete A;
 
  return sum;
}

int Graph::depth_first_search(int u,int v,int dcount,int* low,int* pre) const
{
  int w,output = 0,dc = dcount + 1;
  std::set<int>::const_iterator it;

  pre[v] = dc;
  low[v] = pre[v];
  for(it=neighbours[v].begin(); it!=neighbours[v].end(); ++it) {
    w = *it;
    if (pre[w] == -1) {
      output += depth_first_search(v,w,dc,low,pre);
      low[v] = std::min(low[v],low[w]); 
      if (low[w] == pre[w]) output++;
    }
    else if (w != u) {
      low[v] = std::min(low[v],pre[w]);
    }
  }
  return output;
}

int Graph::bridge_count() const
{
  int i,bcount = 0;
  int low[nvertex],pre[nvertex];

  for(i=0; i<nvertex; ++i) {
    low[i] = -1;
    pre[i] = -1;
  }

  for(i=0; i<nvertex; ++i) {
    if (pre[i] == -1) bcount += depth_first_search(i,i,0,low,pre);
  }
  return bcount;
}

double Graph::cyclicity() const
{
  double output = double(nedge - bridge_count())/double(nedge);
  return output;
}

void Graph::genus(int* chi) const
{
  // Method to compute a collection of topological indices associated with the graph
  // This includes the maximum, minimum and average number of valences, the number of
  // edges, whether the graph is connected, and the maximum and minimum genus. If this
  // is a serial/SMP simulation, we can also easily compute the connectivity and the
  // minimum omega.
  assert(connected());
  chi[0] = int(double(nedge)/6.0 - 0.5*double(nvertex-2));
  chi[1] = int(double((nvertex-3)*(nvertex-4))/12.0);
}

double Graph::connectivity() const 
{
  int i,in1,nn,mm;
  double output = 0.0;
  std::set<int>::const_iterator it;

  // The connectivity is just the sum
  // c = \sum_{edges} 1/sqrt{|v_i|*|v_j|}
  for(i=0; i<nvertex; ++i) {
    for(it=neighbours[i].begin(); it!=neighbours[i].end(); ++it) {
      in1 = *it;
      if (in1 > i) {
        nn = neighbours[i].size();
        mm = neighbours[in1].size();
        output += 1.0/std::sqrt(double(nn*mm));
      }
    }
  }
  return output;
}

int Graph::omega() const 
{
  // The minimum omega is the sum
  // omega = \sum_{i=1}^N 1/(N-|v_i|)
  double sum = 0.0;
  for(int i=0; i<nvertex; ++i) {
    sum += 1.0/double(nvertex - neighbours[i].size());
  }
  return int(sum);
}

void Graph::compute_laplacian(Matrix<double>* laplacian) const
{
  int i;
  std::set<int>::const_iterator it;
  const double m_one = -1.0;

  laplacian->initialize(nvertex,nvertex);

  for(i=0; i<nvertex; ++i) {
    laplacian->set(i,i,double(neighbours[i].size()));
    for(it=neighbours[i].begin(); it!=neighbours[i].end(); ++it) {
      laplacian->set(i,*it,m_one);
    }
  }
}

int Graph::chromatic_number() const
{
  // This method determines the chromatic number of the graph
  int i,j,in1,v,chi;
  bool changed,flag;
  std::set<int>::const_iterator it;
  const int mdegree = max_degree();
  // It makes no sense to try to compute a single chromatic number for 
  // a disconnected graph
  if (!connected()) return 0;

  // The first array contains a histogram of the vertices, ranked by 
  // valence, and the second contains an actual list of the vertices
  int* npopulation = new int[1+mdegree];
  int** population = new int*[1+mdegree];
  // The colour associated with each vertex
  int* colour = new int[nvertex];
  int min_colour;

  // Initialize arrays and finish allocating the memory
  for(i=0; i<=mdegree; ++i) {
    npopulation[i] = 0;
    population[i] = new int[nvertex];
  }
  for(i=0; i<nvertex; ++i) {
    colour[i] = -1;
  }

  // Now we sort the vertices by valence
  for(i=0; i<nvertex; ++i) {
    v = neighbours[i].size();
    population[v][npopulation[v]] = i;
    npopulation[v] += 1;
  }
  // And next move to colouring the vertices, beginning with those 
  // with the greatest valence
  for(i=mdegree; i>0; --i) {
    // Loop through each vertex with this valence
    for(j=0; j<npopulation[i]; ++j) {
      // Try to give it the first colour
      min_colour = 0;
      v = population[i][j];
      // Now loop through all this vertex's neighbours and see which, if 
      // any, colours are already being used. If a colour is in use, then 
      // retry again, from scratch
      flag = true;
      do {
        changed = false;
        for(it=neighbours[v].begin(); it!=neighbours[v].end(); ++it) {
          in1 = *it;
          // This colour is in use, so try the next highest colour and break 
          // out of this loop to retest all the vertex links again, and see 
          // if it doesn't conflict with the neighbours
          if (colour[in1] == min_colour) {
            min_colour++;
            changed = true;
            break;
          }
        }
        if (!changed) flag = false;
      } while(flag);
      // Colour successfully assigned, so now store it in the colour array
      colour[v] = min_colour;
    }
  }

  // The chromatic number is the maximum colour assigned, plus one
  chi = 0;
  for(i=0; i<nvertex; ++i) {
    j = colour[i];
    if (j > chi) chi = j;
  }
  chi += 1;

  // Restore the memory we allocated
  for(i=0; i<=mdegree; ++i) {
    delete[] population[i];
  }
  delete[] population;
  delete[] npopulation;
  delete[] colour;

  return chi;
}

double Graph::compute_energy() const
{
  // To compute the topological energy, first determine if the graph is 
  // connected
  if (!connected()) return 1000000.0;
  int chi[2];
  genus(chi);
  return (double(chi[0])/double(chi[1]));
}

int Graph::minimize_topology(int nsteps,double temperature,std::vector<double>& energy_history)
{
  // The method carries nsteps annealing steps at the requested temperature, storing the 
  // energy at each iteration in the array E_history and returning the number of annealing 
  // steps that were accepted.
  // Note that the only topology changes consist of adding or deleting edges among the 
  // vertices - no attempt is made to alter the number of vertices.
  int l,n,m,acceptance = 0;
  double alpha,boltzmann,E_old,current_energy;
  bool added;

  // Compute the initial geometric energy and store it in E_old
  energy_history.clear();
  E_old = compute_energy();
  energy_history.push_back(E_old);

  // Begin the annealing steps - we can use a for loop here because 
  // a geometric mutation can never disconnect the graph, obviously
  for(l=0; l<nsteps; ++l) {
    // Mutate the topology...
    n = RND.irandom(nvertex);
    do {
      m = RND.irandom(nvertex);
      if (n != m) break;
    } while(true);
    if (neighbours[n].count(m) == 1) {
      // The edge exists, so delete it...
      neighbours[n].erase(m);
      neighbours[m].erase(n);
      nedge--;
      added = false;
    }
    else {
      // This edge doesn't exist, so add it...
      add_edge(n,m);
      added = true;
    }
    // Compute the energy
    current_energy = compute_energy();
    // If the new energy is less than the old energy, accept the step
    if (current_energy < E_old) {
      E_old = current_energy;
      acceptance++;
    }
    else {
      // Otherwise, use the Boltzmann criterion, with CPU #0 doing the 
      // work if this is a distributed/parallel simulation 
      alpha = RND.drandom();
      boltzmann = std::exp((E_old - current_energy)/temperature);
      if (boltzmann > alpha) {
        // In this case, accept the spatial mutation nonetheless, so let the 
        // other CPUs know this, and then increment acceptance and reset the 
        // value of E_old
        E_old = current_energy;
        acceptance++;
      }
      else {
        // This one is rejected, so let the other CPUs know this, and rollback the 
        // mutation
        if (added) {
          neighbours[n].erase(m);
          neighbours[m].erase(n);
          nedge--;
        }
        else {
          add_edge(n,m);
        }
      }
    }
    // Put the current value of E_old in the array E_history
    energy_history.push_back(E_old);
  }
  // Return the number of annealing steps in which the mutation was 
  // accepted
  return acceptance;
}
