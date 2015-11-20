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

#include "graph.h"

using namespace SYNARMOSMA;

extern Random RND;

Graph::Graph() : Schema()
{
  // The empty graph, no vertices or edges...
  nvertex = 0;
}

Graph::Graph(int n,bool complete) : Schema(n)
{
  // If complete = false, just build a graph on n vertices and exit
  if (!complete) return;
  // The complete graph on n vertices...
  int i,j,nc = 0;
  std::set<int> vx;

  for(i=0; i<n; ++i) {
    for(j=1+i; j<n; ++j) {
      vx.insert(i);
      vx.insert(j);
      index_table[vx] = nc;
      edges.push_back(Edge(i,j));
      neighbours[i].insert(j);
      neighbours[j].insert(i);
      vx.clear();
      nc++;
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
        if (rho < double(sum)/double(2*size())) break;
      }
      add_edge(i,k);
    }
  }
}

Graph::Graph(int n,double p) : Schema(n)
{
  // A constructor that builds a graph with n vertices and p percent of 
  // the number of edges in the complete graph on n vertices, chosen randomly.
  int i,j,nc = 0;
  double alpha;
  std::set<int> vx;

  for(i=0; i<n; ++i) {
    for(j=1+i; j<n; ++j) {
      alpha = RND.drandom();
      if (alpha > p) continue;
      vx.insert(i);
      vx.insert(j);
      index_table[vx] = nc;
      edges.push_back(Edge(i,j));
      neighbours[i].insert(j);
      neighbours[j].insert(i);
      vx.clear();
      nc++;
    }
  }
}

Graph::Graph(const Graph& source)
{
  nvertex = source.nvertex;
  neighbours = source.neighbours;
  edges = source.edges;
}

Graph& Graph::operator =(const Graph& source) 
{
  if (this == &source) return *this;
  nvertex = source.nvertex;
  neighbours = source.neighbours;
  edges = source.edges;
  return *this;
}

Graph::~Graph()
{

}

bool Graph::consistent() const
{
  if (!Schema::consistent()) return false;
  int i,j,vx[2],nedge = (signed) edges.size();
  std::set<int>::const_iterator it;
  for(i=0; i<nedge; ++i) {
    if (edges[i].nodes.size() != 2) return false;
    j = 0;
    for(it=edges[i].nodes.begin(); it!=edges[i].nodes.end(); ++it) {
      vx[j] = *it; ++j;
    }
    if (vx[0] < 0 || vx[0] >= nvertex) return false;
    if (vx[1] < 0 || vx[1] >= nvertex) return false;
    if (!edges[i].active) continue;
    if (neighbours[vx[0]].count(vx[1]) == 0) return false;
    if (neighbours[vx[1]].count(vx[0]) == 0) return false;
  }
  return true;
}

void Graph::serialize(std::ofstream& s) const
{
  int i,j,nedge = (signed) edges.size();
  std::set<int>::const_iterator it;

  s.write((char*)(&nvertex),sizeof(int));
  for(i=0; i<nvertex; ++i) {
    j = (signed) neighbours[i].size();
    s.write((char*)(&j),sizeof(int));
    for(it=neighbours[i].begin(); it!=neighbours[i].end(); ++it) {
      j = *it;
      s.write((char*)(&j),sizeof(int));
    }
  }
  s.write((char*)(&nedge),sizeof(int));
  for(i=0; i<nedge; ++i) {
    edges[i].serialize(s);
  }
}

void Graph::deserialize(std::ifstream& s)
{
  int i,j,k,n;
  Edge q;
  std::set<int> S;

  clear();

  s.read((char*)(&nvertex),sizeof(int));
  for(i=0; i<nvertex; ++i) {
    s.read((char*)(&n),sizeof(int));
    for(j=0; j<n; ++j) {
      s.read((char*)(&k),sizeof(int));
      S.insert(k);
    }
    neighbours.push_back(S);
    S.clear();
  }
  s.read((char*)(&n),sizeof(int));
  for(i=0; i<n; ++i) {
    q.deserialize(s);
    edges.push_back(q);
    index_table[q.nodes] = i;
  }
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
    G->edges = edges;
    G->neighbours = neighbours;
    return;
  }

  int i;
  bool found;

  G->nvertex = nvertex; 
  G->edges = edges;
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
        G->drop_vertex(i);
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
  if (size() > (3*nvertex - 6)) return false;
  // Now the hard case where we need to do some work
  // to get the answer...
  bool output = true;

  return output;
}

bool Graph::drop_vertex(int v)
{
  assert(v >= 0 && v < nvertex);
  int i,j,nedges = (signed) edges.size();
  std::vector<int> deletion;
  std::set<int> S = neighbours[v];
  std::set<int>::const_iterator it;

  for(it=S.begin(); it!=S.end(); ++it) {
    neighbours[*it].erase(v);
  }
  //nedge -= S.size();
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
  // In this case we will need to delete the edges containing v, 
  // since the 
  index_table.clear();
  for(i=0; i<nedges; ++i) {
    if (edges[i].nodes.count(v) > 0) {
      deletion.push_back(i);
      continue;
    }
    S.clear();
    for(it=edges[i].nodes.begin(); it!=edges[i].nodes.end(); ++it) {
      j = *it;
      if (j > v) {
        S.insert(j-1);
      }
      else {
        S.insert(j);
      }
    }
    edges[i].nodes = S;
  }
  j = (signed) deletion.size();
  for(i=j-1; i>=0; --i) {
    edges.erase(edges.begin() + deletion[i]);
  }
  for(i=0; i<(signed) edges.size(); ++i) {
    index_table[edges[i].nodes] = i;
  }
  return true;
}

bool Graph::add_edge(int v,int u)
{
  if (u == v) return false;
  if (!connected(v,u)) {
    neighbours[v].insert(u);
    neighbours[u].insert(v);
    std::set<int> vx; vx.insert(v); vx.insert(u);
    hash_map::const_iterator qt = index_table.find(vx);
    if (qt == index_table.end()) {
      index_table[vx] = (signed) edges.size();
      edges.push_back(Edge(v,u));
    }
    else {
      assert(!edges[qt->second].active);
      edges[qt->second].active = true;
    }
    return true;
  }
  return false;
}

bool Graph::drop_edge(int v,int u)
{
  if (u == v) return false;
  if (!connected(v,u)) return false;

  neighbours[u].erase(v);
  neighbours[v].erase(u);
  std::set<int> vx; vx.insert(u); vx.insert(v);
  hash_map::const_iterator qt = index_table.find(vx);
  edges[qt->second].active = false;

  return true;
}

void Graph::clear()
{
  nvertex = 0;
  neighbours.clear();
  edges.clear();
  index_table.clear();
}

int Graph::make_complete()
{
  // A method that outfits the graph with the complete topology, so that 
  // every vertex is connected to every other. It returns the number of edges 
  // that had to be added...
  int i,j,sum = 0;

  for(i=0; i<nvertex; ++i) {
    for(j=1+i; j<nvertex; ++j) {
      if (add_edge(i,j)) sum++;
    }
  }
  return sum;
}

bool Graph::stellar_deletion(int v)
{
  // This only works if this vertex is of degree three...
  if (neighbours[v].size() != 3) return false;
  int n,i = 0,vx[3];
  std::set<int> S = neighbours[v];
  std::set<int>::const_iterator it;
  drop_vertex(v);
  for(it=S.begin(); it!=S.end(); ++it) {
    n = *it;
    if (n > v) n--;
    vx[i] = n;
    ++i;
  }
  add_edge(vx[0],vx[1]);
  add_edge(vx[0],vx[2]);
  add_edge(vx[1],vx[2]);
  return true;
}

bool Graph::stellar_addition(int v)
{
  // This only works if this vertex is part of a 3-cycle...
  int n,m,vx[3];
  std::vector<boost::tuple<int,int,int> > cycles;
  std::set<int> S = neighbours[v];
  std::set<int>::const_iterator it,jt;

  for(it=S.begin(); it!=S.end(); ++it) {
    n = *it;
    for(jt=S.begin(); jt!=S.end(); ++jt) {
      m = *jt;
      if (m == n) continue;
      if (neighbours[n].count(m) > 0) {
        cycles.push_back(boost::make_tuple(v,n,m));
        break;
      }
    }
  }
  if (cycles.empty()) return false;
  m = RND.irandom(cycles.size());
  vx[0] = cycles[m].get<0>();
  vx[1] = cycles[m].get<1>();
  vx[2] = cycles[m].get<2>();
  drop_edge(vx[0],vx[1]);
  drop_edge(vx[0],vx[2]);
  drop_edge(vx[1],vx[2]);
  n = add_vertex();
  add_edge(vx[0],n);
  add_edge(vx[1],n);
  add_edge(vx[2],n);
  return true;
}

bool Graph::fusion(int v,int u)
{
  if (u == v) return false;
  if (u == -1) {
    if (neighbours[v].empty()) return false;
    u = RND.irandom(neighbours[v]);
  }
  if (!drop_vertex(u)) return false;
  int i;
  std::set<int> S = neighbours[u];
  std::set<int>::const_iterator it;
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
  return drop_edge(v,u);
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
  const double NE = double(size());
  const double NV = double(nvertex);
  
  Graph wcopy(*this);

  if (site) {
    // Site percolation - we remove vertices and their associated edges...
    do {
      n = RND.irandom(wcopy.nvertex);
      if (!wcopy.drop_vertex(n)) continue;
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
    output = double(wcopy.size())/NE;
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
  int i,j,p,q,alpha,output = 1 + size(),done[nvertex],parent[nvertex],dist[nvertex];
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
  double ginv = (output == (1+size())) ? 0.0 : 1.0/double(output);
  return ginv;
}

int Graph::cyclomatic_number() const
{
  return (size() - nvertex + 1);
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
  int i,n,output = 1 + size();
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

  output = 2.0*double(size())/double(nvertex*(nvertex-1));
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
  sum -= double(size()); 

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
  const int nedge = size();
  double output = double(nedge - bridge_count())/double(nedge);
  return output;
}

int Graph::genus(std::vector<int>& gamma) const
{
  // A method to compute the genus of a connected graph, assumed simple as well, or to 
  // calculate at least the minimum and maximum genus of this graph.  
  gamma.clear();
  assert(connected());
  if (planar()) return 0;
  int nedge = size();
  int ll = std::ceil(double(nedge)/6.0 - 0.5*double(nvertex) + 1.0);
  int ul = std::floor(0.5*double(nedge - nvertex + 1));
  // If the lower and upper limits are identical, then we have the genus of the graph
  if (ll == ul) return ll;
  gamma.push_back(ll);
  gamma.push_back(ul);
  return -1;
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

void Graph::compute_laplacian(Matrix<double>* L) const
{
  int i;
  std::set<int>::const_iterator it;
  const double m_one = -1.0;

  L->initialize(nvertex,nvertex);

  for(i=0; i<nvertex; ++i) {
    L->set(i,i,double(neighbours[i].size()));
    for(it=neighbours[i].begin(); it!=neighbours[i].end(); ++it) {
      L->set(i,*it,m_one);
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
  std::vector<int> gamma;
  int g = genus(gamma);
  if (g == -1) {
    return double(gamma[0])/(1.0 + double(gamma[1]));
  }
  else {
    return g;
  }
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
    if (connected(n,m)) {
      // The edge exists, so delete it...
      drop_edge(n,m);
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
          drop_edge(n,m);
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
