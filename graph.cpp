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

Graph::~Graph()
{

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

bool Graph::add_edge(int v1,int v2)
{
  if (v1 == v2) return false;
  std::set<int>::const_iterator it = std::find(neighbours[v1].begin(),neighbours[v1].end(),v2);
  if (it == neighbours[v1].end()) {
    neighbours[v1].insert(v2);
    neighbours[v2].insert(v1);
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
      for(it=neighbours[j].begin(); it!=neighbours[j].end(); it++) {
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
  const double nv = double(nvertex);
  int counter[1+max];

  for(i=0; i<=max; ++i) {
    counter[i] = 0;
  }

  for(i=0; i<nvertex; ++i) {
    counter[neighbours[i].size()] += 1;
  }

  histogram.clear();
  histogram.push_back(double(counter[1])/nv);
  if (logarithmic) {
    // Use logarithmic binning, so intervals of size {1,2,4,8,16...}
    int lbound = 1,ubound = 2,sum,its = 1;
    double alpha;
    do {
      lbound *= 2;
      ubound *= 2;
      if (lbound > max) break;
      if (ubound > (1+max)) ubound = 1+max;
      sum = 0;
      for(i=lbound; i<ubound; ++i) {
        sum += counter[i];
      }
      alpha = double(sum)/double(ubound - lbound);
      histogram.push_back(alpha/nv);
      its++;
    } while(true);
  }
  else {
    // Use a uniform bin width of unity...
    for(i=2; i<=max; ++i) {
      histogram.push_back(double(counter[i])/nv);
    }
  }
}

double Graph::percolation(bool site) const 
{
  assert(connected());
  int i,n,m,nc;
  double output;
  std::vector<int> csize,components;
  const double NE = double(nedge);
  const double NV = double(nvertex);
  
  Graph wcopy(*this);

  if (site) {
    // Site percolation - we remove vertices and their associated edges...
    int j;
    std::set<int> S;
    std::set<int>::const_iterator it;

    do {
      n = RND.irandom(wcopy.nvertex);
      S = wcopy.neighbours[n];
      for(it=S.begin(); it!=S.end(); ++it) {
        wcopy.neighbours[*it].erase(n);
      }
      wcopy.nedge -= S.size();
      for(i=0; i<wcopy.nvertex; ++i) {
        if (i == n) continue;
        S.clear();
        for(it=wcopy.neighbours[i].begin(); it!=wcopy.neighbours[i].end(); ++it) {
          j = *it;
          if (j > n) {
            S.insert(j-1);
          }
          else {
            S.insert(j);
          }
        }
        wcopy.neighbours[i] = S;
      }
      wcopy.nvertex--;
      wcopy.neighbours.erase(wcopy.neighbours.begin() + n);
      nc = wcopy.component_analysis(components);
      if (nc == 1) continue;
      // We need to see if the giant component still exists...
      for(i=0; i<nc; ++i) {
        csize.push_back(0);
      }
      for(i=0; i<wcopy.nvertex; ++i) {
        csize[components[i]] += 1;
      }
      std::sort(csize.begin(),csize.end());
      if (double(csize[nc-1])/double(csize[nc-2]) < 0.5) break;
      csize.clear();
    } while(true);
    output = double(wcopy.nvertex)/NV;
  }
  else {
    // Bond percolation - we only remove edges..
    do {
      n = RND.irandom(wcopy.nvertex);
      m = RND.irandom(wcopy.neighbours[n]);
      wcopy.neighbours[n].erase(m);
      wcopy.neighbours[m].erase(n);
      wcopy.nedge--;
      nc = wcopy.component_analysis(components);
      if (nc == 1) continue;
      // We need to see if the giant component still exists...
      for(i=0; i<nc; ++i) {
        csize.push_back(0);
      }
      for(i=0; i<wcopy.nvertex; ++i) {
        csize[components[i]] += 1;
      }
      std::sort(csize.begin(),csize.end());
      if (double(csize[nc-1])/double(csize[nc-2]) < 0.5) break;
      csize.clear();
    } while(true);
    output = double(wcopy.nedge)/NE;
  }

  return output;
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
      for(it=current.begin(); it!=current.end(); it++) {
        p = *it;
        done[p] = 1;
        for(jt=neighbours[p].begin(); jt!=neighbours[p].end(); jt++) {
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

  for(int i=0; i<nvertex*nvertex; ++i) {
    AD[i] = 0.0;
  }

  for(i=0; i<nvertex; ++i) {
    for(it=A->elements[i].begin(); it!=A->elements[i].end(); it++) {
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
    for(it=neighbours[i].begin(); it!=neighbours[i].end(); it++) {
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
