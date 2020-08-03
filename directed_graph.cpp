#include "directed_graph.h"

using namespace SYNARMOSMA;

extern Random RND;

Directed_Graph::Directed_Graph() : Graph() 
{
  
}

Directed_Graph::Directed_Graph(int n) : Graph(n)
{
  int i,j,nc = 0;
  double alpha;
  std::set<int> vx;
  Relation d; 

  for(i=0; i<n; ++i) {
    for(j=1+i; j<n; ++j) {
      vx.insert(i);
      vx.insert(j);
      index_table[vx] = nc;
      d = Relation::disparate;
      alpha = RND.drandom();
      if (alpha < 0.15) {
        d = Relation::before;
        number_directed++;
      }
      else if (alpha < 0.3) {
        d = Relation::after;
        number_directed++;
      }      
      edges.push_back(Edge(i,j,0.0,d));
      neighbours[i].insert(j);
      neighbours[j].insert(i);
      vx.clear();
      nc++;
    }
  }
}

Directed_Graph::Directed_Graph(int n,double p) : Graph(n)
{
  if (p < -std::numeric_limits<double>::epsilon()) throw std::invalid_argument("The edge probability must be greater than zero!");
  if (p > 1.0) throw std::invalid_argument("The edge probability must not be greater than one!");
  int i,j,nc = 0;
  double alpha;
  std::set<int> vx;
  auto d = Relation::disparate;

  number_directed = 0;

  for(i=0; i<n; ++i) {
    for(j=1+i; j<n; ++j) {
      alpha = RND.drandom();
      if (alpha > p) continue;
      vx.insert(i);
      vx.insert(j);
      index_table[vx] = nc;
      d = Relation::disparate;
      alpha = RND.drandom();
      if (alpha < 0.15) {
        d = Relation::before;
        number_directed++;
      }
      else if (alpha < 0.3) {
        d = Relation::after;
        number_directed++;
      }      
      edges.push_back(Edge(i,j,0.0,d));
      neighbours[i].insert(j);
      neighbours[j].insert(i);
      vx.clear();
      nc++;
    }
  }
}

Directed_Graph::~Directed_Graph()
{

}

Directed_Graph::Directed_Graph(const Directed_Graph& source)
{
  nvertex = source.nvertex;
  neighbours = source.neighbours;
  edges = source.edges;
  index_table = source.index_table;
  number_directed = source.number_directed;
}

Directed_Graph& Directed_Graph::operator =(const Directed_Graph& source) 
{
  if (this == &source) return *this;

  nvertex = source.nvertex;
  neighbours = source.neighbours;
  edges = source.edges;
  index_table = source.index_table;
  number_directed = source.number_directed;

  return *this;
}

int Directed_Graph::serialize(std::ofstream& s) const
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
  s.write((char*)(&number_directed),sizeof(int)); count += sizeof(int);

  return count;
}

int Directed_Graph::deserialize(std::ifstream& s) 
{
  int i,j,k,n,vx[2],count = 0;
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
    q.get_vertices(vx);
    S.clear();
    S.insert(vx[0]); S.insert(vx[1]);
    index_table[S] = i;
  }
  s.read((char*)(&number_directed),sizeof(int)); count += sizeof(int);

  return count;

}

void Directed_Graph::write2disk(const std::string& filename) const
{
  int i,j;
  Relation rho;
  std::set<int> S;
  std::set<int>::const_iterator it;
  hash_map::const_iterator qt;

  std::ofstream s(filename,std::ios::trunc);

  s << "digraph G {" << std::endl;
  // First all the elements in the poset...
  for(i=0; i<nvertex; ++i) {
    s << "  \"" << 1+i << "\";" << std::endl;
  }
  // Now the directed edges induced by the poset's ordering...
  for(i=0; i<nvertex; ++i) {
    for(it=neighbours[i].begin(); it!=neighbours[i].end(); ++it) {
      j = *it;
      S.clear();
      S.insert(i); S.insert(j);
      qt = index_table.find(S);
      rho = edges[qt->second].get_direction(i,j);
      if (rho == Relation::before) {
        s << "  \"" << 1+i << "\" -> \"" << 1+j << "\";" << std::endl;
      }
      else if (rho == Relation::after) {
        s << "  \"" << 1+j << "\" -> \"" << 1+i << "\";" << std::endl;
      }
      else {
        s << "  \"" << 1+i << "\" -- \"" << 1+j << "\";" << std::endl;
      }
    }
  }
  s << "}" << std::endl;
  s.close();
}

int Directed_Graph::maximum_parents() const
{
  // Find the maximum number of parents over the entire set of vertices of this directed 
  // graph...
  int i,j,nparents,max_parents = -1;
  std::set<int> S;
  std::set<int>::const_iterator jt;
  hash_map::const_iterator qt;

  for(i=0; i<nvertex; ++i) {
    nparents = 0;
    for(jt=neighbours[i].begin(); jt!=neighbours[i].end(); ++jt) {
      j = *jt;
      S.clear();
      S.insert(i); S.insert(j);
      qt = index_table.find(S);
      // Does this edge point from j to i?
      if (edges[qt->second].get_direction(j,i) == Relation::before) nparents++;
    }
    if (nparents > max_parents) max_parents = nparents;    
  }

  return max_parents;
}

bool Directed_Graph::singly_connected() const
{
  // This method tests whether or not for every pair of vertices i and j in the graph, there is 
  // no more than one path leading from i to j.
  int i,j,k,u,v,npath;
  bool visited[nvertex];
  std::set<int> S,current,next;
  std::set<int>::const_iterator it,jt;
  hash_map::const_iterator qt;

  for(i=0; i<nvertex; ++i) {
    for(j=0; j<nvertex; ++j) {
      if (i == j) continue;
      // Check if it is possible to get from i to j 
      for(k=0; k<nvertex; ++k) {
        visited[k] = false;
      }
      npath = 0;
      current.clear();
      next.clear();
      current.insert(i);
      visited[i] = true;
      do {
        for(it=current.begin(); it!=current.end(); ++it) {
          u = *it;
          for(jt=neighbours[u].begin(); jt!=neighbours[u].end(); ++jt) {
            v = *jt;
            if (visited[v]) continue;
            S.clear();
            S.insert(u); S.insert(v);
            qt = index_table.find(S);
            if (edges[qt->second].get_direction(u,v) == Relation::before) next.insert(v);
          }
        }
        if (next.empty()) break;
        current.clear();
        for(it=next.begin(); it!=next.end(); ++it) {
          if (*it == j) {
            npath++;
          }
          else {
            visited[*it] = true;
            current.insert(*it);
          }
        }
        next.clear();
      } while(true);
      if (npath > 1) return false;
    }
  }
  return true;
}

bool Directed_Graph::connected() const
{
  int i,n = -1,p,m = 0;
  std::vector<int> ubiquity;
  std::set<int> S,change,nchange;
  std::set<int>::const_iterator it,jt;
  hash_map::const_iterator qt;

  for(i=0; i<nvertex; ++i) {
    p = (signed) neighbours[i].size();
    if (p > m) {
      m = p;
      n = i;
    }
    ubiquity.push_back(0);
  }
  if (n == -1) return false;
  ubiquity[n] = 1;
  change.insert(n);

  do {
    for(it=change.begin(); it!=change.end(); ++it) {
      i = *it;
      for(jt=neighbours[i].begin(); jt!=neighbours[i].end(); ++jt) {
        n = *jt;
        if (ubiquity[n] == 1) continue;
        S.clear();
        S.insert(i); S.insert(n);
        qt = index_table.find(S);
        if (edges[qt->second].get_direction(i,n) == Relation::before) nchange.insert(n);
      }
    }
    if (nchange.empty()) break;
    for(it=nchange.begin(); it!=nchange.end(); ++it) {
      ubiquity[*it] = 1;
    }
    change = nchange;
    nchange.clear();
  } while(true);
  for(i=0; i<nvertex; ++i) {
    if (ubiquity[i] == 0) return false;
  }

  return true;
}

int Directed_Graph::distance(int u,int v) const
{
  // A method to calculate the topological distance from 
  // u to v; the method returns -1 if it is impossible to 
  // get from u to v.
  if (u == v) return 0;

  int i,j,delta = -1,its = 1;
  bool visited[nvertex];
  std::set<int> S,current,next;
  std::set<int>::const_iterator it,jt;
  hash_map::const_iterator qt;

  for(i=0; i<nvertex; ++i) {
    visited[i] = false;
  }
  current.insert(u);
  visited[u] = true;
  do {
    for(it=current.begin(); it!=current.end(); ++it) {
      i = *it;
      for(jt=neighbours[i].begin(); jt!=neighbours[i].end(); ++jt) {
        j = *jt;
        if (visited[j]) continue;
        S.clear();
        S.insert(i); S.insert(j);
        qt = index_table.find(S);
        if (edges[qt->second].get_direction(i,j) == Relation::before) next.insert(j);
      }
    }
    if (next.empty()) break;
    if (next.count(v) > 0) {
      delta = its;
      break;
    }
    for(it=next.begin(); it!=next.end(); ++it) {
      visited[*it] = true;
    }
    current = next;
    next.clear();
    its++;
  } while(true);

  return delta;
}

void Directed_Graph::compute_distances(pair_index& output) const
{
  int i,j,k,delta;
  std::pair<int,int> pr;
  std::vector<int> distances;
  std::set<int> S;
  std::set<int>::const_iterator it;
  hash_map::const_iterator qt;

  for(i=0; i<nvertex; ++i) {
    for(j=0; j<nvertex; ++j) {
      distances.push_back(std::numeric_limits<int>::max()); 
    }
  }
  for(i=0; i<nvertex; ++i) {
    distances[nvertex*i+i] = 0;
    for(it=neighbours[i].begin(); it!=neighbours[i].end(); ++it) {
      j = *it;
      S.clear();
      S.insert(i); S.insert(j);
      qt = index_table.find(S);
      if (edges[qt->second].get_direction(i,j) == Relation::before) {      
        distances[nvertex*i+j] = 1;
      }
    }
  }

  for(k=0; k<nvertex; ++k) {
    for(i=0; i<nvertex; ++i) {
      if (distances[nvertex*i+k] == std::numeric_limits<int>::max()) continue;
      for(j=0; j<nvertex; ++j) {
        if (distances[nvertex*k+j] == std::numeric_limits<int>::max()) continue;
        delta = distances[nvertex*i+k] + distances[nvertex*k+j];
        if (delta < distances[nvertex*i+j]) distances[nvertex*i+j] = delta;
      }
    }
  }
  output.clear();
  for(i=0; i<nvertex; ++i) {
    for(j=0; j<nvertex; ++j) {
      if (i == j) continue;
      if (distances[nvertex*i+j] == std::numeric_limits<int>::max()) distances[nvertex*i+j] = -1;
      pr.first = i; pr.second = j;
      output[pr] = distances[nvertex*i+j];
      pr.first = j; pr.second = i;
      output[pr] = distances[nvertex*j+i];
    }
  }
}

bool Directed_Graph::add_edge(int u,int v)
{
  return add_edge(u,v,Relation::disparate,0.0);
}

bool Directed_Graph::add_edge(int u,int v,double kappa)
{
  return add_edge(u,v,Relation::disparate,kappa);
}

bool Directed_Graph::add_edge(int u,int v,Relation d,double ell)
{
  bool output = Graph::add_edge(u,v,ell);

  std::set<int> S;
  S.insert(u); S.insert(v);
  hash_map::const_iterator qt = index_table.find(S);

  if (d == Relation::disparate) {
    edges[qt->second].set_direction(d);
  }
  else {
    if (u < v) {
      edges[qt->second].set_direction(d);
    }
    else {
      Relation rho = (d == Relation::before) ? Relation::after : Relation::before;
      edges[qt->second].set_direction(rho);
    }
  }
  return output;
}

bool Directed_Graph::mutate_edge(int u,int v)
{
  if (u == v) return false;
  std::set<int> S;
  S.insert(u); S.insert(v);
  hash_map::const_iterator qt = index_table.find(S);
  if (qt == index_table.end()) return false;
  int n = qt->second;
  double alpha = RND.drandom();
  auto d = edges[n].get_direction();
  if (d == Relation::disparate) {
    d = (alpha < 0.5) ? Relation::before : Relation::after;
  }
  else {
    if (d == Relation::before) {
      d = (alpha < 0.25) ? Relation::after : Relation::disparate;
    }
    else {
      d = (alpha < 0.25) ? Relation::before : Relation::disparate;
    }
  }
  edges[n].set_direction(d);

  return true;
}

bool Directed_Graph::directed_cycle(const std::vector<int>& path,int base,int length) const
{
  if (base == path[0] && path.size() > 2) return true;
  if ((signed) path.size() == length) return false;
  int v;
  std::set<int> S;
  std::set<int>::const_iterator it;
  hash_map::const_iterator qt;

  for(it=neighbours[base].begin(); it!=neighbours[base].end(); ++it) {
    v = *it;
    S.clear();
    S.insert(base); S.insert(v);
    qt = index_table.find(S);
    if (edges[qt->second].get_direction(base,v) == Relation::before) {
      std::vector<int> npath = path;
      npath.push_back(v);
      if (directed_cycle(npath,v,length)) return true;
    }
  }
  return false;
}

bool Directed_Graph::acyclic() const
{
  std::vector<int> path;

  for(int i=0; i<nvertex; ++i) {
    path.clear();
    path.push_back(i);
    if (directed_cycle(path,i,nvertex)) return false;
  }
  return true;
}

double Directed_Graph::compute_flow(int source,int sink) 
{
  int i;
  bool valid = false;
  std::set<int> S;
  std::set<int>::const_iterator it;
  hash_map::const_iterator qt;
  const int ne = size();

  for(i=0; i<ne; ++i) {
    if (edges[i].get_direction() == Relation::disparate) edges[i].set_capacity(0.0);
  }
  // Next we need to verify that there is at least one outgoing edge with capacity > 0 
  // for the source and at least one incoming edge with capacity > 0 for the sink
  for(it=neighbours[source].begin(); it!=neighbours[source].end(); ++it) {
    i = *it;
    S.clear();
    S.insert(source);
    S.insert(i);
    qt = index_table.find(S);
    if (edges[qt->second].get_capacity() < std::numeric_limits<double>::epsilon()) continue;
    if (edges[qt->second].get_direction(source,i) == Relation::before) valid = true;
    if (valid) break;
  }
  if (!valid) {
    for(i=0; i<ne; ++i) {
      edges[i].set_flow(0.0);
    }
    return 0.0;
  }
  // Now the sink
  valid = false;
  for(it=neighbours[sink].begin(); it!=neighbours[sink].end(); ++it) {
    i = *it;
    S.clear();
    S.insert(sink);
    S.insert(i);
    qt = index_table.find(S);
    if (edges[qt->second].get_capacity() < std::numeric_limits<double>::epsilon()) continue;
    if (edges[qt->second].get_direction(sink,i) == Relation::after) valid = true;
    if (valid) break;
  }
  if (!valid) {
    for(i=0; i<ne; ++i) {
      edges[i].set_flow(0.0);
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
    if (edges[i].get_direction() == Relation::before) {
      pr.first = vx[0]; pr.second = vx[1];
      rgraph[pr] = int(10000.0*edges[i].get_capacity());
    }
    else if (edges[i].get_direction() == Relation::after) {
      pr.first = vx[1]; pr.second = vx[0];
      rgraph[pr] = int(10000.0*edges[i].get_capacity());
    }
  }

  int max_flow = network_flow(rgraph,source,sink,nvertex);
  for(i=0; i<ne; ++i) {
    edges[i].set_flow(0.0);
    edges[i].get_vertices(vx);
    if (edges[i].get_direction() == Relation::before) {
      pr.first = vx[0]; pr.second = vx[1]; 
      edges[i].set_flow(edges[i].get_capacity() - double(rgraph[pr])/10000.0);
    }
    else if (edges[i].get_direction() == Relation::after) {
      pr.first = vx[1]; pr.second = vx[0];
      edges[i].set_flow(edges[i].get_capacity() - double(rgraph[pr])/10000.0);
    }
  }
  return double(max_flow)/10000.0;
}

unsigned int Directed_Graph::compute_sinks(std::set<int>& output) const
{
  int i,j,test;
  std::set<int> S;
  std::set<int>::const_iterator it;
  hash_map::const_iterator qt;

  output.clear();
  // A sink is a vertex all of whose edges are incoming or undirected...
  for(i=0; i<nvertex; ++i) {
    if (neighbours[i].empty()) continue;
    test = 0;
    for(it=neighbours[i].begin(); it!=neighbours[i].end(); ++it) {
      j = *it;
      S.clear();
      S.insert(i); S.insert(j);
      qt = index_table.find(S);
      if (edges[qt->second].get_direction(i,j) == Relation::before) {
        // If this vertex has an outgoing edge, it isn't a sink...
        test = -1;
        break;
      }
      else if (edges[qt->second].get_direction(i,j) == Relation::after) {
        test = 1;
      }
    }
    // To be a sink, the vertex must have at least one incoming edge and no outgoing edges
    if (test == 1) output.insert(i);
  }
  return output.size();
}

unsigned int Directed_Graph::compute_sources(std::set<int>& output) const
{
  int i,j,test;
  std::set<int> S;
  std::set<int>::const_iterator it;
  hash_map::const_iterator qt;

  output.clear();
  // A source is a vertex all of whose edges are outgoing or undirected...
  for(i=0; i<nvertex; ++i) {
    if (neighbours[i].empty()) continue;
    test = 0;
    for(it=neighbours[i].begin(); it!=neighbours[i].end(); ++it) {
      j = *it;
      S.clear();
      S.insert(i); S.insert(j);
      qt = index_table.find(S);
      if (edges[qt->second].get_direction(i,j) == Relation::after) {
        // If this vertex has an incoming edge, it isn't a source...
        test = -1;
        break;
      }
      else if (edges[qt->second].get_direction(i,j) == Relation::before) {
        test = 1;
      }
    }
    // To be a source, the vertex must have at least one outgoing edge and no incoming edges
    if (test == 1) output.insert(i);
  }
  return output.size();
}
