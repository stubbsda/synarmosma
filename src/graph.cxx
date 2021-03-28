#include "graph.h"

using namespace SYNARMOSMA;

Graph::Graph()
{
  // The empty graph, no vertices or edges...
}

Graph::Graph(int n) : Schema(n)
{
  // A graph with n vertices but no edges...
  if (n < 1) throw std::invalid_argument("The graph order must be greater than zero!");
}

Graph::Graph(const std::string& name)
{
  // An assortment of "named" graphs...
  std::string uname = boost::to_upper_copy(name);
  if (uname == "DÜRER") {
    for(int i=0; i<12; ++i) {
      add_vertex();
    }
    add_edge(0,2);
    add_edge(0,4);
    add_edge(0,6);

    add_edge(1,3);
    add_edge(1,5);
    add_edge(1,7);

    add_edge(2,4);
    add_edge(2,8);

    add_edge(3,5);
    add_edge(3,9);

    add_edge(4,10);

    add_edge(5,11);

    add_edge(6,7);
    add_edge(6,11);

    add_edge(7,8);

    add_edge(8,9);

    add_edge(9,10);

    add_edge(10,11);
  }
  else if (uname == "GOLOMB") {
    for(int i=0; i<10; ++i) {
      add_vertex();
    }
    add_edge(0,1);
    add_edge(0,5);
    add_edge(0,6);

    add_edge(1,2);
    add_edge(1,6);
    add_edge(1,7);

    add_edge(2,3);
    add_edge(2,6);

    add_edge(3,4);
    add_edge(3,6);
    add_edge(3,8);

    add_edge(4,5);
    add_edge(4,6);

    add_edge(5,6);
    add_edge(5,9);

    add_edge(7,8);
    add_edge(7,9);
    add_edge(8,9);
  }
  else if (uname == "HERSCHEL") {
    for(int i=0; i<11; ++i) {
      add_vertex();
    }
    add_edge(0,1);
    add_edge(0,7);
    add_edge(0,8);
    add_edge(0,10);

    add_edge(1,2);
    add_edge(1,9);

    add_edge(2,3);
    add_edge(2,10);

    add_edge(3,4);
    add_edge(3,9);

    add_edge(4,5);
    add_edge(4,8);
    add_edge(4,10);

    add_edge(5,6);
    add_edge(5,9);

    add_edge(6,7);
    add_edge(6,8);

    add_edge(7,9);
  }
  else if (uname == "PETERSEN") {
    for(int i=0; i<10; ++i) {
      add_vertex();
    }
    add_edge(0,1);
    add_edge(0,4);
    add_edge(0,5);

    add_edge(1,2);
    add_edge(1,6);
   
    add_edge(2,3);
    add_edge(2,7);

    add_edge(3,4);
    add_edge(3,8);
  
    add_edge(4,9);

    add_edge(5,7);
    add_edge(5,8);

    add_edge(6,8);
    add_edge(6,9);
        
    add_edge(7,9);
  }
  else if (uname == "WAGNER") {
    for(int i=0; i<8; ++i) {
      add_vertex();
    }
    for(int i=0; i<7; ++i) {
      add_edge(i,i+1);
    }
    add_edge(0,7);

    for(int i=0; i<4; ++i) {
      add_edge(i,i+4);
    }
  }
  else {
    throw std::invalid_argument("Unrecognized name in graph constructor!");
  }
}

Graph::Graph(int n,const std::string& type) : Schema(n)
{
  if (n < 1) throw std::invalid_argument("The graph order must be greater than zero!");

  std::string utype = boost::to_upper_copy(type);
  if (utype == "COMPLETE") {
    // The complete graph on n vertices...
    for(int i=0; i<n; ++i) {
      for(int j=1+i; j<n; ++j) {
        add_edge(i,j);
      }
    }
  }
  else if (utype == "CHAIN") {
    // A minimally connected graph with n - 1 edges
    for(int i=0; i<n-1; ++i) {
      add_edge(i,i+1);
    }
  }
  else if (utype == "CYCLIC") {
    // Much like the chain model except with a cyclic topology, 
    // thus a final edge connecting the end of the chain to its 
    // beginning
    for(int i=0; i<n-1; ++i) {
      add_edge(i,i+1);
    }
    // The final edge that makes it cyclic
    add_edge(0,n-1); 
  }
  else if (utype == "CONNECTED") {
    // Create a random connected graph on n vertices - we keep adding 
    // random edges until the graph is connected
    if (n < 2) throw std::invalid_argument("The order of a connected graph must be greater than one!");

    int u,v;
    Random RND;

    do {
      do {
        u = RND.irandom(n);
        v = RND.irandom(n);
        if (u != v) break;
      } while(true);
      if (!add_edge(u,v)) continue;
      if (connected()) break;
    } while(true);
  }
  else {
    throw std::invalid_argument("Unrecognized type in graph constructor!");
  }
}

Graph::Graph(int n,int c) : Schema(n)
{
  if (n < 1) throw std::invalid_argument("The graph order must be greater than zero!");
  if (c < 1) throw std::invalid_argument("The minimum degree in the graph constructor must be greater than zero!");

  // Construct a scale-free graph with n vertices where the minimum degree is c
  int i,j,k,sum;
  double rho;
  Random RND;

  for(i=0; i<c+1; ++i) {
    for(j=(signed) neighbours[i].size(); j<c; ++j) {
      do {
        k = RND.irandom(c+1);
        if (k == i) continue;
        if (neighbours[i].count(k) == 0) break;
      } while(true);
      add_edge(i,k,0.0);
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
      add_edge(i,k,0.0);
    }
  }
}

Graph::Graph(int n,double p) : Schema(n)
{
  if (n < 1) throw std::invalid_argument("The graph order must be greater than zero!");
  if (p < std::numeric_limits<double>::epsilon() || p > 1.0) throw std::invalid_argument("The edge percentage must lie between 0 and 1!");

  // A constructor that builds a graph with n vertices and p percent of 
  // the number of edges in the complete graph on n vertices, chosen randomly.
  int i,j;
  double alpha;
  Random RND;

  for(i=0; i<n; ++i) {
    for(j=1+i; j<n; ++j) {
      alpha = RND.drandom();
      if (alpha > p) continue;
      add_edge(i,j,0.0);
    }
  }
}

Graph::Graph(int n,const std::vector<int>& btype)
{
  build_lattice(n,btype);
}

Graph::Graph(const Graph& source)
{
  nvertex = source.nvertex;
  neighbours = source.neighbours;
  edges = source.edges;
  index_table = source.index_table;
}

Graph& Graph::operator =(const Graph& source) 
{
  if (this == &source) return *this;

  nvertex = source.nvertex;
  neighbours = source.neighbours;
  edges = source.edges;
  index_table = source.index_table;

  return *this;
}

Graph::~Graph()
{

}

void Graph::initialize(const Graph& G)
{
  nvertex = G.nvertex;
  neighbours = G.neighbours;
  edges = G.edges;
  index_table = G.index_table;  
}

int Graph::build_lattice(int n,const std::vector<int>& btype)
{
  if (n < 2) throw std::invalid_argument("The number of vertices per dimension in Graph::build_lattice must be greater than one!");
  int d = (signed) btype.size();
  if (d < 2) throw std::invalid_argument("The number of dimensions in Graph::build_lattice must be greater than one!");
  for(int i=0; i<d; ++i) {
    if (btype[i] != 0 && btype[i] != 1) throw std::invalid_argument("The graph's boundary topology in Graph::build_lattice must be either linear or toroidal!");
  }

  // This constructor builds a graph with the topology of a rectangular or Cartesian lattice, of 
  // dimensionality equal to the length of the vector btype and with n vertices per dimension. The 
  // btype vector determines the nature of the boundary for each dimension: 0 means it is linear, 
  // 1 means toroidal.
  int i,j,k,p,q,in1,base;
  std::set<int> S;
  std::vector<int> vx,index;

  clear();

  nvertex = int(ipow(n,d));
  for(i=0; i<nvertex; ++i) {
    neighbours.push_back(S);
  }

  for(i=0; i<d; ++i) {
    index.push_back(0);
  }

  for(i=0; i<nvertex; ++i) {
    q = i;
    for(j=0; j<d; ++j) {
      base = ipow(n,d-1-j);
      p = q/base;
      index[j] = p;
      q = q - p*base;
    }
    for(j=0; j<d; ++j) {
      vx = index;
      if (index[j] > 0) {
        vx[j] -= 1;
        in1 = 0;
        for(k=0; k<d; ++k) {
          in1 += ipow(n,d-1-k)*vx[k];
        }
        if (i < in1) add_edge(i,in1); 
      }
      else {
        if (btype[j] == 1) {
          vx[j] = n-1;
          in1 = 0;
          for(k=0; k<d; ++k) {
            in1 += ipow(n,d-1-k)*vx[k];
          }
          if (i < in1) add_edge(i,in1); 
        }
      }
      vx = index;
      if (index[j] < n-1) {
        vx[j] += 1;
        in1 = 0;
        for(k=0; k<d; ++k) {
          in1 += ipow(n,d-1-k)*vx[k];
        }
        if (i < in1) add_edge(i,in1); 
      }
    }
  }
  return size();
}

bool Graph::consistent() const
{
  if (!Schema::consistent()) return false;

  int i,vx[2];
  std::set<int> S;
  hash_map::const_iterator qt;
  const int nedge = (signed) edges.size();
 
  for(i=0; i<nedge; ++i) {
    edges[i].get_vertices(vx);
    // Edge is a self-loop...
    if (vx[0] == vx[1]) return false;
    // Edge vertices don't exist...
    if (vx[0] < 0 || vx[0] >= nvertex) return false;
    if (vx[1] < 0 || vx[1] >= nvertex) return false;
    // Edge is inconsistent with neighbour table...
    if (neighbours[vx[0]].count(vx[1]) == 0) return false;
    if (neighbours[vx[1]].count(vx[0]) == 0) return false;
    S.clear();
    S.insert(vx[0]); S.insert(vx[1]); 
    qt = index_table.find(S);
    // Unable to find this edge in the edge index...
    if (qt == index_table.end()) return false;
    // Wrong index value...
    if (qt->second != i) return false;
  }
  return true;
}

int Graph::serialize(std::ofstream& s) const
{
  int i,j,count = 0,nedge = (signed) edges.size();
  std::set<int>::const_iterator it;

  s.write((char*)(&nvertex),sizeof(int)); count += sizeof(int);
  for(i=0; i<nvertex; ++i) {
    j = (signed) neighbours[i].size();
    s.write((char*)(&j),sizeof(int)); count += sizeof(int);
    for(it=neighbours[i].begin(); it!=neighbours[i].end(); ++it) {
      j = *it;
      s.write((char*)(&j),sizeof(int)); count += sizeof(int);
    }
  }
  s.write((char*)(&nedge),sizeof(int)); count += sizeof(int);
  for(i=0; i<nedge; ++i) {
    count += edges[i].serialize(s);
  }
  return count;
}

int Graph::deserialize(std::ifstream& s)
{
  int i,j,k,n,count = 0;
  Edge q;
  std::set<int> S;

  clear();

  s.read((char*)(&nvertex),sizeof(int)); count += sizeof(int);
  for(i=0; i<nvertex; ++i) {
    s.read((char*)(&n),sizeof(int)); count += sizeof(int);
    for(j=0; j<n; ++j) {
      s.read((char*)(&k),sizeof(int)); count += sizeof(int);
      S.insert(k);
    }
    neighbours.push_back(S);
    S.clear();
  }
  s.read((char*)(&n),sizeof(int)); count += sizeof(int);
  for(i=0; i<n; ++i) {
    count += q.deserialize(s);
    edges.push_back(q);
    S.clear();
    S.insert(q.low); S.insert(q.high);
    index_table[S] = i;
  }
  return count;
}

int Graph::compute_surface(const std::set<int>& vx,std::set<int>& surface) const
{
  // This method returns the number of nodes which lie on the boundary of 
  // the set of nodes contained in the method's argument. A node is on the 
  // boundary if it has at least one edge which joins it to a vertex which 
  // is not a member of the set "vx". This method will also compute the set 
  // of surface edges, i.e. the set of edges which connect the nodes of vx 
  // to the rest of the graph.
  int i,j;
  std::set<int> S,bnodes;
  hash_map::const_iterator qt;
  std::set<int>::const_iterator it,jt;

  surface.clear();
 
  for(it=vx.begin(); it!=vx.end(); ++it) {
    i = *it;
    for(jt=neighbours[i].begin(); jt!=neighbours[i].end(); ++jt) {
      j = *jt;
      if (vx.count(j) > 0) continue;
      bnodes.insert(i);
      S.clear();
      S.insert(i); S.insert(j);
      qt = index_table.find(S);
      surface.insert(qt->second);
    }
  }
  return (signed) bnodes.size();
}

int Graph::compute_hamiltonian_path(bool cycle,int nattempts,std::vector<int>& path,int start) const
{
  if (!connected()) {
    path.clear();
    return -1;
  }
  int i,u,v,l,n,svertex,winner;
  bool visited[nvertex],solution;
  std::vector<std::pair<int,int> > candidates;
  std::set<int> next;
  std::set<int>::const_iterator it,jt;
  Random RND;

  for(l=0; l<nattempts; ++l) {
    svertex = (start == -1) ? RND.irandom(nvertex) : start;
    path.clear();
    path.push_back(svertex);
    for(i=0; i<nvertex; ++i) {
      visited[i] = false;
    }
    visited[svertex] = true;
    do {
      v = path.back();
      for(it=neighbours[v].begin(); it!=neighbours[v].end(); ++it) {
        u = *it;
        if (visited[u]) continue;
        n = 0;
        for(jt=neighbours[u].begin(); jt!=neighbours[u].end(); ++jt) {
          if (!visited[*jt]) n++;
        } 
        candidates.push_back(std::pair<int,int>(u,n));
      }
      if (candidates.empty()) break;
      n = (signed) candidates.size();
      std::sort(candidates.begin(),candidates.end(),pair_predicate_int);
      u = candidates[0].second;
      for(i=1; i<n; ++i) {
        if (candidates[i].second == u) next.insert(candidates[i].first);
      }
      if (next.empty()) {
        winner = candidates[0].first;
      }
      else {
        next.insert(candidates[0].first);
        winner = RND.irandom(next);
        next.clear();
      }
      path.push_back(winner);
      visited[winner] = true;
      if ((signed) path.size() == nvertex) break;
      candidates.clear();
    } while(true);
    // Check if this is a valid path, i.e. one which covers the whole graph...
    if ((signed) path.size() < nvertex) continue;
    solution = true;
    for(i=0; i<nvertex; ++i) {
      if (!visited[i]) {
        solution = false;
        break;
      }
    }
    if (!solution) continue;
    // So this is a Hamiltonian path, now we need to check if it is a cycle, if necessary...
    if (!cycle) break; 
    u = path.back();
    if (neighbours[u].count(svertex) > 0) break;
    solution = false;
  }

  if (!solution) path.clear();

  return 1+l;
}

void Graph::write2disk(const std::string& filename,const std::vector<std::string>& names) const
{
  int i,j;
  std::set<int>::const_iterator it;

  std::ofstream s(filename,std::ios::trunc);

  s << "graph G {" << std::endl;
  if (names.empty()) {
    // First all the vertices...
    for(i=0; i<nvertex; ++i) {
      s << "  \"" << 1+i << "\";" << std::endl;
    }
    // Now the edges...
    for(i=0; i<nvertex; ++i) {
      for(it=neighbours[i].begin(); it!=neighbours[i].end(); ++it) {
        j = *it;
        if (i > j) continue;
        s << "  \"" << 1+i << "\" -- \"" << 1+j << "\";" << std::endl;
      }
    }
  }
  else {
    // First all the vertices...
    for(i=0; i<nvertex; ++i) {
      s << "  \"" << names[i] << "\";" << std::endl;
    }
    // Now the edges...
    for(i=0; i<nvertex; ++i) {
      for(it=neighbours[i].begin(); it!=neighbours[i].end(); ++it) {
        j = *it;
        if (i > j) continue;
        s << "  \"" << names[i] << "\" -- \"" << names[j] << "\";" << std::endl;
      }
    }
  }
  s << "}" << std::endl;
  s.close();
}

void Graph::core(Graph* G,int k) const
{
  // A method to compute the k-core of the current graph
  if (k < 1) throw std::invalid_argument("For the k-core of a graph, k must be greater than zero!");

  G->clear();

  // The k-core in this case is null
  if (k > max_degree()) return;

  G->initialize(*this);

  // The 1-core and the k-core when k < min_degree(G) is just the graph itself...
  if (k == 1 || k < min_degree()) return;

  // The less trivial cases...
  int i,n;

  do {
    if (k > G->max_degree()) {
      // The k-core will necessarily be null in this case, so we can exit
      G->clear();
      break;
    }   
    n = G->nvertex;
    for(i=0; i<n; ++i) {
      if (k > G->degree(i)) {
        G->drop_vertex(i);
        break;
      }
    }
    if (n != G->nvertex) break;
  } while(true);
}

void Graph::vertex_centrality(std::vector<double>& output,double tolerance,bool katz) const
{
  // A method to compute the centrality of the vertices in the 
  // graph - we first need to compute the adjacency matrix and its 
  // largest eigenvalue, then we use a recursive process in order to 
  // compute the vector of centrality values. This method can use 
  // either the standard formula or that for Katz centrality.
  if (nvertex < 2) throw std::runtime_error("The vertex centrality can only be computed if the graph order is greater than one!");
  int i;
  double sum;
  std::vector<double> x,xnew;
  const double beta = 1.0;

  adjacency_eigenvalues(x);

  const double alpha = (katz) ? 0.9/x[nvertex-1] : 0.5/x[nvertex-1];

  x.clear();
  for(i=0; i<nvertex; ++i) {
    x.push_back(1.0);
    xnew.push_back(0.0);
  }

  Binary_Matrix* A = new Binary_Matrix;
  compute_adjacency_matrix(A);

  do {
    A->multiply(x,xnew);
    for(i=0; i<nvertex; ++i) {
      xnew[i] *= alpha;
      xnew[i] += beta;
    }
    sum = 0.0;
    for(i=0; i<nvertex; ++i) {
      sum += (xnew[i] - x[i])*(xnew[i] - x[i]);
    }
    if (std::sqrt(sum) < tolerance) break;
    x = xnew;
  } while(true);

  delete A;

  output = xnew;
}

int Graph::DFS_bridge(int u,int v,int dcount,int* low,int* pre) const
{
  int w,output = 0,dc = dcount + 1;
  std::set<int>::const_iterator it;

  pre[v] = dc;
  low[v] = pre[v];
  for(it=neighbours[v].begin(); it!=neighbours[v].end(); ++it) {
    w = *it;
    if (pre[w] == -1) {
      output += DFS_bridge(v,w,dc,low,pre);
      low[v] = std::min(low[v],low[w]); 
      if (low[w] == pre[w]) output++;
    }
    else if (w != u) {
      low[v] = std::min(low[v],pre[w]);
    }
  }

  return output;
}

int Graph::DFS_cycle(int u,int v,std::vector<int>& path,bool* visited) const
{
  int w,out = 0;
  std::set<int>::const_iterator it;

  visited[v] = true;
  path.push_back(v);
  for(it=neighbours[v].begin(); it!=neighbours[v].end(); ++it) {
    w = *it;
    if (!visited[w]) {
      std::vector<int> npath = path;
      out = DFS_cycle(v,w,npath,visited);
    }
    else if (w != u) {
      int n = (signed) path.size();
      for(int i=0; i<n; ++i) {
        if (w == path[i]) return (n-i);
      }
    }
  }
  return out;
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
    if (pre[i] == -1) bcount += DFS_bridge(i,i,0,low,pre);
  }
  return bcount;
}

bool Graph::bipartite() const
{
  // We need to determine if this graph has a cycle of odd length
  int i,length;
  bool visited[nvertex];
  std::vector<int> path;
  
  for(i=0; i<nvertex; ++i) {
    visited[i] = false;
  }

  for(i=0; i<nvertex; ++i) {
    if (!visited[i]) {
      length = DFS_cycle(i,i,path,visited);
      if (length%2 != 0) return false;
    }
    path.clear();
  }
  return true;
}

int Graph::compactness(int v,int n) const
{
  if (v < 0 || v >= nvertex) throw std::invalid_argument("Illegal vertex argument in the Graph::compactness method!");
  if (n < 0) throw std::invalid_argument("Illegal number of steps argument in the Graph::compactness method!");
  
  // The simplest case...
  if (n == 0) return 1;
  int l = 1 + (signed) neighbours[v].size();
  // Check if the number of steps is one or if this graph is complete, in which case we're already done...
  if (n == 1 || l == order()) return l;
  int i,p,q;
  std::set<int> next,current,vx = neighbours[v];
  std::set<int>::const_iterator it,jt;

  current = vx;
  vx.insert(v);

  for(i=2; i<=n; ++i) {
    for(it=current.begin(); it!=current.end(); ++it) {
      p = *it;
      for(jt=neighbours[p].begin(); jt!=neighbours[p].end(); ++jt) {
        q = *jt;
        if (vx.count(q) > 0) continue;
        next.insert(q);
      }
    }
    if (next.empty()) break;
    for(it=next.begin(); it!=next.end(); ++it) {
      vx.insert(*it);
    }
    current = next;
    next.clear();
  }
  return (signed) vx.size();
}

int Graph::eccentricity(int v) const
{
  if (v < 0 || v >= nvertex) throw std::invalid_argument("The Graph::eccentricity argument must be a vertex of the graph!");

  int i,delta,output = 0;
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i,delta) reduction(max:output)
#endif
  for(i=0; i<v; ++i) {
    delta = distance(i,v);
    if (delta > output) output = delta;
  }
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i,delta) reduction(max:output)
#endif
  for(i=v+1; i<nvertex; ++i) {
    delta = distance(i,v);
    if (delta > output) output = delta;
  }
  return output;  
}

int Graph::radius() const
{
  if (!connected()) return -1;
  int i,t,output = size();

  for(i=0; i<nvertex; ++i) {
    t = eccentricity(i);
    if (t < output) output = t;
  }
  return output;
}

int Graph::diameter() const
{
  if (!connected()) return -1;
  int i,t,output = -1;

  for(i=0; i<nvertex; ++i) {
    t = eccentricity(i);
    if (t > output) output = t;
  }
  return output;
}

bool Graph::planar() const
{
  if (!connected()) throw std::runtime_error("To compute a graph's planarity it must be connected!");

  // Small graphs are always planar...
  if (nvertex <= 2) return true;
  // The more edges, the less likely it's planar...
  if (nvertex > 2 && size() > (3*nvertex - 6)) return false;
  // If every vertex has degree greater than five, it can't be planar...
  if (min_degree() > 5) return false;
  int i,vx[2],ne = size();
  boost::adjacency_list<boost::vecS,boost::vecS,boost::undirectedS,boost::property<boost::vertex_index_t,int> > G(nvertex);
  for(i=0; i<ne; ++i) {
    edges[i].get_vertices(vx);
    boost::add_edge(vx[0],vx[1],G);
  }
  return boyer_myrvold_planarity_test(G);
}

bool Graph::drop_vertex(int v)
{
  if (v < 0 || v >= nvertex) throw std::invalid_argument("The Graph::drop_vertex argument must be a vertex of the graph!");

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
  // We now need to delete the edges containing v
  index_table.clear();
  for(i=0; i<nedges; ++i) {
    if (edges[i].low == v || edges[i].high == v) {
      deletion.push_back(i);
      continue;
    }
    if (edges[i].low > v) {
      edges[i].low -= 1;
    }
    if (edges[i].high > v) {
      edges[i].high -= 1;
    }
  }
  j = (signed) deletion.size();
  for(i=j-1; i>=0; --i) {
    edges.erase(edges.begin() + deletion[i]);
  }
  for(i=0; i<(signed) edges.size(); ++i) {
    S.clear();
    S.insert(edges[i].low); S.insert(edges[i].high);
    index_table[S] = i;
  }
  return true;
}

bool Graph::add_edge(int v,int u)
{
  return add_edge(v,u,0.0);
}

bool Graph::add_edge(int v,int u,double kappa)
{
  if (u == v) return false;
  std::set<int> vx; 
  vx.insert(v); vx.insert(u);
  if (!connected(v,u)) {
    neighbours[v].insert(u);
    neighbours[u].insert(v);
    index_table[vx] = (signed) edges.size();
    edges.push_back(Edge(v,u,kappa));
    return true;
  }
  hash_map::const_iterator qt = index_table.find(vx);
  edges[qt->second].length = kappa;
  return false;
}

bool Graph::drop_edge(int v,int u)
{
  if (u == v) return false;
  if (!connected(v,u)) return false;

  std::set<int> vx; 
  hash_map::const_iterator qt;

  vx.insert(u); vx.insert(v);
  qt = index_table.find(vx);
  // This edge doesn't exist!
  if (qt == index_table.end()) return false;
  neighbours[u].erase(v);
  neighbours[v].erase(u);
  edges.erase(edges.begin() + qt->second);
  index_table.clear();
  for(int i=0; i<(signed) edges.size(); ++i) {
    vx.clear();
    vx.insert(edges[i].low); vx.insert(edges[i].high);
    index_table[vx] = i;
  }
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
  if (!drop_vertex(v)) return false;
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
  std::vector<std::tuple<int,int,int> > cycles;
  std::set<int> S = neighbours[v];
  std::set<int>::const_iterator it,jt;
  Random RND;

  for(it=S.begin(); it!=S.end(); ++it) {
    n = *it;
    for(jt=S.begin(); jt!=S.end(); ++jt) {
      m = *jt;
      if (m == n) continue;
      if (neighbours[n].count(m) > 0) {
        cycles.push_back(std::make_tuple(v,n,m));
        break;
      }
    }
  }
  if (cycles.empty()) return false;
  m = RND.irandom(cycles.size());
  vx[0] = std::get<0>(cycles[m]);
  vx[1] = std::get<1>(cycles[m]);
  vx[2] = std::get<2>(cycles[m]);
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
    Random RND;
    u = RND.irandom(neighbours[v]);
  }
  std::set<int> S = neighbours[u];
  if (!drop_vertex(u)) return false;
  int i;
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
  if (u == -1) {
    Random RND;
    u = RND.irandom(neighbours[v]);
  }
  return drop_edge(v,u);
}

bool Graph::foliation_m(int v,int u)
{
  if (u == -1) {
    Random RND;
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

double Graph::compute_flow(int source,int sink) 
{
  int i;
  bool valid = false;
  std::set<int> S;
  std::set<int>::const_iterator it;
  hash_map::const_iterator qt;
  const int ne = size();

  // Next we need to verify that there is at least one edge with capacity > 0 
  // for the source and at least one edge with capacity > 0 for the sink
  for(it=neighbours[source].begin(); it!=neighbours[source].end(); ++it) {
    S.clear();
    S.insert(source);
    S.insert(*it);
    qt = index_table.find(S);
    if (edges[qt->second].capacity > std::numeric_limits<double>::epsilon()) {
      valid = true;
      break;
    }
  }
  if (!valid) {
    for(i=0; i<ne; ++i) {
      edges[i].flow = 0.0;
    }
    return 0.0;
  }
  // Now the sink
  valid = false;
  for(it=neighbours[sink].begin(); it!=neighbours[sink].end(); ++it) {
    S.clear();
    S.insert(sink);
    S.insert(*it);
    qt = index_table.find(S);
    if (edges[qt->second].capacity > std::numeric_limits<double>::epsilon()) {
      valid = true;
      break;
    }
  }
  if (!valid) {
    for(i=0; i<ne; ++i) {
      edges[i].flow = 0.0;
    }
    return 0.0;
  }
  // So we can finally proceed to computing a non-trivial flow on this directed graph
  int vx[2];
  std::pair<int,int> pr;
  pair_index rgraph;
  // Create a residual graph and fill the residual graph with
  // given capacities in the original graph as residual capacities
  // in residual graph
  for(i=0; i<ne; ++i) {
    edges[i].get_vertices(vx);
    pr.first = vx[0]; pr.second = vx[1];
    rgraph[pr] = int(10000.0*edges[i].capacity);
    pr.first = vx[1]; pr.second = vx[0];
    rgraph[pr] = int(10000.0*edges[i].capacity);
  }
  int max_flow = network_flow(rgraph,source,sink,nvertex);
  for(i=0; i<ne; ++i) {
    edges[i].flow = 0.0;
    edges[i].get_vertices(vx);
    pr.first = vx[0]; pr.second = vx[1];
    if (rgraph[pr] > 0) {      
      edges[i].flow = edges[i].capacity - double(rgraph[pr])/10000.0;
    }
    pr.first = vx[1]; pr.second = vx[0];
    if (rgraph[pr] > 0) {
      edges[i].flow = -(edges[i].capacity - double(rgraph[pr])/10000.0);
    }
  }
  return double(max_flow)/10000.0;
}

void Graph::complement(Graph* G) const
{
  int i,j;
  std::set<int> S;
  hash_map::const_iterator qt;

  G->clear();
  for(i=0; i<nvertex; ++i) {
    G->add_vertex();
  }
  for(i=0; i<nvertex; ++i) {
    for(j=1+i; j<nvertex; ++j) {
      S.clear();
      S.insert(i); S.insert(j);
      qt = index_table.find(S);
      if (qt == index_table.end()) G->add_edge(i,j);
    }
  }
}

void Graph::defoliate(const Pseudograph* parent,std::vector<Monomial<int> >& tutte) const
{
  std::vector<int> cvector;
  int nb = parent->get_candidates(cvector);
  Random RND;

  if (cvector.empty()) {
    int nl = parent->self_loops();
    Monomial<int> term;
    term.coefficient = 1;
    if (nb > 0) term.exponents.push_back(std::pair<unsigned int,unsigned int>(0,nb));
    if (nl > 0) term.exponents.push_back(std::pair<unsigned int,unsigned int>(1,nl));
    tutte.push_back(term);
    return;
  }
  int cd = (signed) cvector.size()/2;
  int e = RND.irandom(cd);
  int u = cvector[2*e];
  int v = cvector[2*e+1];

  Pseudograph* c1 = new Pseudograph(*parent);
  parent->remove(u,v,c1);
  defoliate(c1,tutte);
  delete c1;

  Pseudograph* c2 = new Pseudograph(*parent);
  parent->contract(u,v,c2);
  defoliate(c2,tutte);
  delete c2;
}

void Graph::tutte_polynomial(std::vector<Monomial<int> >& output) const
{
  int i,j,cf,nt;
  unsigned int k,p,q,p1,p2,base;
  Monomial<int> term;
  std::set<int>::const_iterator it;
  std::vector<Monomial<int> > tutte;
  Pseudograph* G = new Pseudograph(nvertex);

  for(i=0; i<nvertex; ++i) {
    for(it=neighbours[i].begin(); it!=neighbours[i].end(); ++it) {
      j = *it;
      if (i > j) continue;
      G->add_edge(i,j,Relation::disparate);
    }
  }
  defoliate(G,tutte);
  delete G; 
  nt = (signed) tutte.size();
  bool occupied[nt];
  for(i=0; i<nt; ++i) {
    occupied[i] = false;
  }
  output.clear();
  for(i=0; i<nt; ++i) {
    if (occupied[i]) continue;
    term = tutte[i];
    cf = term.coefficient;
    k = tutte[i].exponents.size();
    if (k == 1) {
      base = tutte[i].exponents[0].first;
      if (base == 0) {
        p1 = tutte[i].exponents[0].second; p2 = 0;
      }
      else {
        p1 = 0; p2 = tutte[i].exponents[0].second;
      }
    }
    else {
      p1 = tutte[i].exponents[0].second;
      p2 = tutte[i].exponents[1].second;
    }  
    for(j=1+i; j<nt; ++j) {
      k = tutte[j].exponents.size();
      if (k == 1) {
        base = tutte[j].exponents[0].first;
        if (base == 0) {
          p = tutte[j].exponents[0].second; q = 0;
        }
        else {
          p = 0; q = tutte[j].exponents[0].second;
        }
      }
      else {
        p = tutte[j].exponents[0].second;
        q = tutte[j].exponents[1].second;
      }
      if (p1 == p && p2 == q) {
        cf += 1;
        occupied[j] = true;
      }  
    }
    occupied[i] = true;
    term.coefficient = cf;
    output.push_back(term);
  }
}

double Graph::clustering_coefficient(int v) const
{
  // This method calculates the percentage of distinct pairs (u,w) 
  // of neighbours of v which are also connected directly
  if (neighbours[v].size() < 2) return 0.0;
  unsigned int i,j,n,np,nf = 0;
  std::vector<int> vx;
  std::set<int>::const_iterator it;

  for(it=neighbours[v].begin(); it!=neighbours[v].end(); ++it) {
    vx.push_back(*it);
  }
  n = vx.size();
  for(i=0; i<n; ++i) {
    for(j=1+i; j<n; ++j) {
      if (neighbours[vx[i]].count(vx[j]) > 0 && neighbours[vx[j]].count(vx[i]) > 0) nf++;
    }
  }
  np = n*(n - 1)/2;
  return double(nf)/double(np);
}

double Graph::clustering_coefficient() const
{
  int i;
  double sum = 0.0;

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i) reduction(+:sum)
#endif
  for(i=0; i<nvertex; ++i) {
    sum += clustering_coefficient(i);
  }
  return sum/double(nvertex);
}

double Graph::mean_path_length() const
{
  int i,j;
  double sum = 0.0,npair = double(nvertex*(nvertex-1)/2);

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i,j) reduction(+:sum) schedule(dynamic,1)
#endif
  for(i=0; i<nvertex; ++i) {
    for(j=1+i; j<nvertex; ++j) {
      sum += double(distance(i,j));
    }
  }
  sum = sum/npair;
  return sum;
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

double Graph::entropy(int n) const
{
  if (n < 0) throw std::invalid_argument("The argument of Graph::entropy must be positive!");
  if (size() == 0) std::runtime_error("The entropy of a graph with no edges is undefined!");
  if (n == 1) return std::log(size());
  if (n == 0) n = nvertex - 1;

  int i,j;
  double sum = 0.0;
  std::set<int>::const_iterator it;
  Matrix<double> A(nvertex),B(nvertex);

  for(i=0; i<nvertex; ++i) {
    for(it=neighbours[i].begin(); it!=neighbours[i].end(); ++it) {
      j = *it;
      if (i < j) {
        A.set(i,j,1.0);
        A.set(j,i,1.0);
      }
    }
  }
  B = A;
  for(i=0; i<n-1; ++i) {
    B = A*B;
  }
  for(i=0; i<nvertex; ++i) {
    for(j=1+i; j<nvertex; ++j) {
      sum += B.get(i,j);
    }
  }

  return std::log(sum);
}

double Graph::return_probability(int base,int length,int ntrials) const
{
  int i,j,next,current,sum = 0;
  unsigned long s = (unsigned) std::time(nullptr);
  bool home;

#ifdef _OPENMP
#pragma omp parallel default(shared) private(i,j,next,current,home) firstprivate(s)
  {
  s *= long(1 + omp_get_thread_num());
#endif
  Random RND;

  RND.set_seed(s);
#pragma omp for reduction(+:sum) schedule(dynamic,1)
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
    if (home) sum++;
  }

#ifdef _OPENMP
  }
#endif

  double rho = double(sum)/double(ntrials);
  return rho;
}

std::pair<double,double> Graph::random_walk(int L,double p,int ntrials) const
{
  int i,j,v;
  double mu,rho = 0.0,sigma = 0.0;
  std::set<int> vx;
  const int nbase = int(p*nvertex);
  double value[nbase];
  Random RND;

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

  std::pair<double,double> output(rho,sigma);
  return output;
}

void Graph::degree_distribution(bool logarithmic,std::vector<double>& histogram) const
{
  if (!connected()) throw std::invalid_argument("It is meaningless to compute the degree distribution for a disconnected graph!");

  int i;
  const int m = max_degree();
  int counter[1+m];

  for(i=0; i<=m; ++i) {
    counter[i] = 0;
  }

  for(i=0; i<nvertex; ++i) {
    counter[degree(i)] += 1;
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
      if (lbound > m) break;
      if (ubound > (1 + m)) ubound = 1 + m;
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
    for(i=2; i<=m; ++i) {
      histogram.push_back(double(counter[i]));
    }
  }
}

double Graph::percolation(bool site) const 
{
  if (!connected()) throw std::invalid_argument("Percolation computations are meaningless for a disconnected graph!");

  int i,n,nc;
  bool success;
  std::vector<int> csize,components;
  Graph* G = new Graph;
  Random RND;
  const double NE = double(size());
  const double NV = double(nvertex);
  
  G->initialize(*this);

  do {
    n = RND.irandom(G->order());
    // Site percolation (we remove vertices and their associated edges) or bond percolation (we only remove edges)...
    success = (site) ? G->drop_vertex(n) : G->foliation_x(n,-1);
    if (!success) continue;
    nc = G->component_analysis(components);
    if (nc == 1) continue;
    // We need to see if the giant component still exists...
    for(i=0; i<nc; ++i) {
      csize.push_back(0);
    }
    for(i=0; i<G->order(); ++i) {
      csize[components[i]] += 1;
    }
    // Sort the component sizes in ascending order
    std::sort(csize.begin(),csize.end());
    // The giant component has vanished if the largest component is less than 
    // twice the size of the next largest component...
    if (double(csize[nc-1])/double(csize[nc-2]) < 2.0) break;
    csize.clear();
  } while(true);

  double output = (site) ? double(G->order())/NV : double(G->size())/NE;

  delete G;

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
  const double denominator = std::sqrt(double(degree(u)*degree(v)));
  compute_adjacency_matrix(A);
  // Now we need to square this binary matrix A...
  for(i=0; i<nvertex; ++i) {
    if (A->get(u,i) && A->get(i,v)) sum += 1;
  }
  delete A;
  return double(sum)/denominator;
}

int Graph::girth() const
{
  if (nvertex < 1) throw std::invalid_argument("The girth cannot be calculated for a graph with no vertices!");

  int length = 1 + size();
#ifdef _OPENMP
#pragma omp parallel default(shared)
  {
#endif
  int i,j,p,q,alpha,done[nvertex],parent[nvertex],dist[nvertex];
  std::set<int> current,next;
  std::set<int>::const_iterator it,jt;

#ifdef _OPENMP
#pragma omp for
#endif
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
#ifdef _OPENMP
#pragma omp critical
            {
#endif
            if (alpha < length) length = alpha;
#ifdef _OPENMP
            }
#endif
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
#ifdef _OPENMP
  }
#endif
  // The graph is acyclic...
  if (length == (1 + size())) return -1;
  return length;
}

double Graph::algebraic_connectivity() const
{
  if (!connected() || nvertex < 2) throw std::invalid_argument("The algebraic connectivity cannot be computed for a disconnected graph or one with an order less than two!");

  int i,j,info,nv = nvertex,nwork = 3*nvertex - 1;
  char jtype = 'N';
  char uplo = 'U';
  double A[nvertex*nvertex],w[nvertex],work[nwork],output = 0.0;  

  for(i=0; i<nvertex*nvertex; ++i) {
    A[i] = 0.0;
  }

  for(i=0; i<nvertex; ++i) {
    A[(nvertex+1)*i] = double(degree(i));
    for(j=1+i; j<nvertex; ++j) {
      if (!connected(i,j)) continue;
      A[nvertex*i+j] = -1.0; 
      A[nvertex*j+i] = -1.0;
    }
  }

  dsyev_(&jtype,&uplo,&nv,A,&nv,w,work,&nwork,&info);
  if (info == 0) output = w[1];
  
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
      A->set(i,j); 
      A->set(j,i); 
    }
  }
}

void Graph::adjacency_eigenvalues(std::vector<double>& output) const
{
  if (!connected() || nvertex < 2) throw std::invalid_argument("The adjacency matrix eigenvalues cannot be computed for a disconnected graph or one with an order less than two!");

  int i,j,info,nv = nvertex,nwork = 3*nvertex - 1;
  char jtype = 'N';
  char uplo = 'U';
  double AD[nvertex*nvertex],w[nvertex],work[nwork];

  output.clear();

  for(i=0; i<nvertex*nvertex; ++i) {
    AD[i] = 0.0;
  }

  for(i=0; i<nvertex; ++i) {
    for(j=1+i; j<nvertex; ++j) {
      if (!connected(i,j)) continue;
      AD[nvertex*i+j] = 1.0; 
      AD[nvertex*j+i] = 1.0;
    }
  }

  // If info is zero after this LAPACK call, then w will contain
  // the eigenvalues of A (which are real since A is symmetric and
  // real) in ascending order, thus w[nv-1] will be the largest
  dsyev_(&jtype,&uplo,&nv,AD,&nv,w,work,&nwork,&info);
  if (info != 0) throw std::runtime_error("Unable to compute the eigenvalues of this graph's adjacency matrix!");

  for(i=0; i<nvertex; ++i) {
    output.push_back(w[i]);
  }
}

double Graph::median_degree() const
{
  if (nvertex < 2) return 0.0;

  const int N = max_degree();
  int i,dcount[N+1],psum,sum = 0;

  for(i=0; i<=N; ++i) {
    dcount[i] = 0;
  }
  for(i=0; i<nvertex; ++i) {
    dcount[neighbours[i].size()] += 1;
  }
  for(i=0; i<=N; ++i) {
    psum = sum;
    sum += dcount[i];
    // A regular graph or one with an exact median, so exit...
    if (dcount[i] == nvertex || 2*sum == nvertex) return double(i);
    // Need to compute the corrected median...
    if (2*sum > nvertex) break;
  }
  double output = double(i - 1) + (0.5*double(nvertex) - double(psum))/double(sum - psum);
  return output;
}

double Graph::entwinement() const
{
  // This method produces a real number between 0 and 1 that measures the
  // degree of "labyrinthicity" of the graph
  if (nvertex < 3) return 0.0;
  // The complete graph has an entwinement of unity, of course...
  if (size() == nvertex*(nvertex-1)/2) return 1.0;

  double output = 0.0;
  std::vector<double> w;

  adjacency_eigenvalues(w);

  // We know that w[nv-1] <= max_degree(), so we divide by this
  // to normalize the output
  if (!w.empty()) {
    double ds = sqrt(double(max_degree()));
    double da = average_degree();
    double m1 = (ds < da) ? da : ds;
    output = (w[nvertex-1] - m1)/(ds*ds - m1);
  }
  return output;
}

double Graph::cyclic_resistance() const
{
  if (!connected()) throw std::invalid_argument("The cyclic resistance is undefined for a disconnected graph!");

  int i,j,k,info,nv = nvertex,nwork = 5*nvertex;
  bool rvector[nvertex];
  double sum = 0.0;
  std::vector<double> delta;
  Matrix<double>* L = new Matrix<double>;
  Matrix<double>* W = new Matrix<double>;
  Binary_Matrix* A = new Binary_Matrix;
  double* C = new double[nvertex*nvertex];
  double* work = new double[nwork];
  int* pivots = new int[nvertex];

  compute_adjacency_matrix(A);
  compute_laplacian(L);
  L->convert(C,'r');
  // Use an LAPACK routine to compute the inverse of the graph Laplacian matrix
  dgetri_(&nv,C,&nv,pivots,work,&nwork,&info);

  W->initialize(nvertex,nvertex);
  for(i=0; i<nvertex; ++i) {
    for(j=i+1; j<nvertex; ++j) {
      sum = 1.0/(C[nvertex*i+i] + C[nvertex*j+j] - (C[nvertex*i+j] + C[nvertex*j+i]));
      W->set(i,j,sum);
      W->set(j,i,sum);
    }
  }
  L->initialize(nvertex,nvertex);
  for(i=0; i<nvertex; ++i) {
    for(j=0; j<nvertex; ++j) {
      sum = 0.0;
      A->get_row(j,rvector);
      for(k=0; k<nvertex; ++k) {
        if (rvector[k]) sum += W->get(k,i);
      }
      L->set(i,j,sum); 
    }
  }
  L->get_diagonal(delta);

  for(i=0; i<nvertex; ++i) {
    sum += delta[i];
  }
  sum -= double(size()); 

  delete[] work;
  delete[] pivots;
  delete[] C;

  delete A;
  delete L;
  delete W;

  return sum;
}

double Graph::cyclicity() const
{
  const int nedge = size();
  double output = (nedge == 0) ? 0.0 : double(nedge - bridge_count())/double(nedge);
  return output;
}

int Graph::genus(std::vector<int>& gamma) const
{
  // A method to compute the genus of a connected graph, assumed simple as well, or to 
  // calculate at least the minimum and maximum genus of this graph.  
  if (!connected()) throw std::invalid_argument("The genus is undefined for a disconnected graph!");

  gamma.clear();
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
  int i,j;
  double output = 0.0;
  std::set<int>::const_iterator it;

  // The connectivity is just the sum
  // c = \sum_{edges} 1/sqrt{|v_i|*|v_j|}
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i,j,it) reduction(+:output)
#endif
  for(i=0; i<nvertex; ++i) {
    for(it=neighbours[i].begin(); it!=neighbours[i].end(); ++it) {
      j = *it;
      if (j > i) output += 1.0/std::sqrt(double(degree(i)*degree(j)));
    }
  }
  return output;
}

int Graph::compute_laplacian(Matrix<double>* L) const
{
  int i,nzero = 0;
  std::set<int>::const_iterator it;

  L->initialize(nvertex,nvertex);

  for(i=0; i<nvertex; ++i) {
    if (neighbours[i].empty()) continue;
    L->set(i,i,double(degree(i))); nzero++;
    for(it=neighbours[i].begin(); it!=neighbours[i].end(); ++it) {
      L->set(i,*it,-1.0); nzero++;
    }
  }
  return nzero;
}

int Graph::compute_deformed_laplacian(std::complex<double> s,Matrix<std::complex<double> >* L) const
{
  int i,nzero = 0;
  std::set<int>::const_iterator it;

  L->initialize(nvertex,nvertex);

  for(i=0; i<nvertex; ++i) {
    if (neighbours[i].empty()) continue;
    L->set(i,i,1.0 + s*s*double(degree(i) - 1)); nzero++;
    for(it=neighbours[i].begin(); it!=neighbours[i].end(); ++it) {
      L->set(i,*it,-s); nzero++;
    }
  }
  return nzero;
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
  if (!connected()) return std::numeric_limits<double>::infinity();
  std::vector<int> gamma;
  int g = genus(gamma);
  double output = (g == -1) ? double(gamma[0])/(1.0 + double(gamma[1])) : double(g);
  return output;
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
  Random RND;

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

namespace SYNARMOSMA {
  Graph operator *(const Graph& G,const Graph& H)
  {
    if (G.order() != H.order()) throw std::invalid_argument("The tensor product operands must have the same order!");
    const int N = G.order(); 
    if (N < 2) throw std::invalid_argument("The tensor product operands must have an order greater than one!");
    const int N2 = N*N;

    int i,j,g1,g2,h1,h2;
    Graph output(N2);

    for(i=0; i<N2; ++i) {
      g1 = i/N;
      h1 = i - N*g1;
      for(j=1+i; j<N2; ++j) {
        g2 = j/N;
        h2 = j - N*g2;
        if (G.connected(g1,g2) && H.connected(h1,h2)) output.add_edge(i,j);
      }
    }
    return output;
  }
}
