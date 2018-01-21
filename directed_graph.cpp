#include "directed_graph.h"

using namespace SYNARMOSMA;

extern Random RND;

Directed_Graph::Directed_Graph() : Graph() 
{
  number_directed = 0;
}

Directed_Graph::Directed_Graph(int n) : Graph(n)
{
  int i,j,d,nc = 0;
  double alpha;
  std::set<int> vx;

  number_directed = 0;

  for(i=0; i<n; ++i) {
    for(j=1+i; j<n; ++j) {
      vx.insert(i);
      vx.insert(j);
      index_table[vx] = nc;
      d = UNDIRECTED;
      alpha = RND.drandom();
      if (alpha < 0.15) {
        d = OUTGOING;
        number_directed++;
      }
      else if (alpha < 0.3) {
        d = INCOMING;
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
  int i,j,d,nc = 0;
  double alpha;
  std::set<int> vx;

  number_directed = 0;

  for(i=0; i<n; ++i) {
    for(j=1+i; j<n; ++j) {
      alpha = RND.drandom();
      if (alpha > p) continue;
      vx.insert(i);
      vx.insert(j);
      index_table[vx] = nc;
      d = UNDIRECTED;
      alpha = RND.drandom();
      if (alpha < 0.15) {
        d = OUTGOING;
        number_directed++;
      }
      else if (alpha < 0.3) {
        d = INCOMING;
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
  s.read((char*)(&number_directed),sizeof(int)); count += sizeof(int);

  return count;

}

void Directed_Graph::compute_directedness() 
{
  // A method to calculate how many of the edges are directed...
  int null = 0;
  std::vector<Edge>::const_iterator it;

  for(it=edges.begin(); it!=edges.end(); ++it) {
    if (it->direction == UNDIRECTED) null++;
  }
  number_directed = size() - null;
}

int Directed_Graph::distance(int u,int v) const
{
  // A method to calculate the topological distance from 
  // u to v; the method returns -1 if it is impossible to 
  // get from u to v.
  if (u == v) return 0;

  int i,j,l = -1,its = 1;
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
        if (edges[qt->second].get_direction(i,j) == OUTGOING) {
          next.insert(j);
        }
      }
    }
    if (next.empty()) break;
    if (next.count(v) > 0) {
      l = its;
      break;
    }
    for(it=next.begin(); it!=next.end(); ++it) {
      visited[*it] = true;
    }
    current = next;
    next.clear();
    its++;
  } while(true);
  return l;
}

void Directed_Graph::compute_distances(edge_hash& output) const
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
      if (edges[qt->second].get_direction(i,j) == OUTGOING) {      
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

bool Directed_Graph::add_edge(int u,int v,int d,double ell)
{
  if (!Graph::add_edge(u,v,ell)) return false;
  if (d != UNDIRECTED) { 
    std::set<int> S;
    S.insert(u); S.insert(v);
    hash_map::const_iterator qt = index_table.find(S);
    edges[qt->second].direction = (u < v) ? d : -d;
  }
  return true;
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
  int d = edges[n].direction;
  if (d == UNDIRECTED) {
    edges[n].direction = (alpha < 0.5) ? OUTGOING : INCOMING;
  }
  else {
    if (d == OUTGOING) {
      edges[n].direction = (alpha < 0.25) ? INCOMING : UNDIRECTED;
    }
    else {
      edges[n].direction = (alpha < 0.25) ? OUTGOING : UNDIRECTED;
    }
  }
  return true;
}

bool Directed_Graph::path_connected(int u,int v) const
{
  // Is it possible to get from u to v following the orientation of the graph edges?
  int i,j,visited[nvertex];
  bool output = false;
  std::set<int> current,next,S;
  std::set<int>::const_iterator it,jt;
  hash_map::const_iterator qt;

  for(i=0; i<nvertex; ++i) {
    visited[i] = 0;
  }

  current.insert(u);
  visited[u] = 1;
  do {
    for(it=current.begin(); it!=current.end(); ++it) {
      i = *it;
      for(jt=neighbours[i].begin(); jt!=neighbours[i].end(); ++jt) {
        j = *jt;
        if (visited[j] == 1) continue;
        S.clear();
        S.insert(i); S.insert(j);
        qt = index_table.find(S);
        if (edges[qt->second].get_direction(i,j) == OUTGOING) next.insert(j);
      }
    }
    if (next.empty()) break;
    if (next.count(v) > 0) {
      output = true;
      break;
    }
    for(it=next.begin(); it!=next.end(); ++it) {
      visited[*it] = 1;
    }
    current = next;
    next.clear();
  } while(true);
  return output;
}

bool Directed_Graph::directed_cycle(const std::vector<int>& path,int base,int length) const
{
  if (base == path[0] && path.size() > 2) return true;
  if ((signed) path.size() == length) return false;
  int v;
  bool out;
  std::set<int> S;
  std::set<int>::const_iterator it;
  hash_map::const_iterator qt;

  for(it=neighbours[base].begin(); it!=neighbours[base].end(); ++it) {
    v = *it;
    S.clear();
    S.insert(base); S.insert(v);
    qt = index_table.find(S);
    if (edges[qt->second].get_direction(base,v) == OUTGOING) {
      std::vector<int> npath = path;
      npath.push_back(v);
      out = directed_cycle(npath,v,length);
      if (out) return true;
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
    if (edges[i].direction == UNDIRECTED) edges[i].capacity = 0.0;
  }
  // Next we need to verify that there is at least one outgoing edge with capacity > 0 
  // for the source and at least one incoming edge with capacity > 0 for the sink
  for(it=neighbours[source].begin(); it!=neighbours[source].end(); ++it) {
    i = *it;
    S.clear();
    S.insert(source);
    S.insert(i);
    qt = index_table.find(S);
    if (edges[qt->second].capacity < std::numeric_limits<double>::epsilon()) continue;
    if (edges[qt->second].get_direction(source,i) == OUTGOING) valid = true;
    if (valid) break;
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
    i = *it;
    S.clear();
    S.insert(sink);
    S.insert(i);
    qt = index_table.find(S);
    if (edges[qt->second].capacity < std::numeric_limits<double>::epsilon()) continue;
    if (edges[qt->second].get_direction(sink,i) == INCOMING) valid = true;
    if (valid) break;
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
  edge_hash rgraph;
  edge_hash::const_iterator qtt;
  // Create a residual graph and fill the residual graph with
  // given capacities in the original graph as residual capacities
  // in residual graph
  for(i=0; i<ne; ++i) {
    edges[i].get_vertices(vx);
    if (edges[i].direction == OUTGOING) {
      pr.first = vx[0]; pr.second = vx[1];
      rgraph[pr] = int(10000.0*edges[i].capacity);
    }
    else if (edges[i].direction == INCOMING) {
      pr.first = vx[1]; pr.second = vx[0];
      rgraph[pr] = int(10000.0*edges[i].capacity);
    }
  }

  int max_flow = network_flow(rgraph,source,sink,nvertex);
  for(i=0; i<ne; ++i) {
    edges[i].flow = 0.0;
    edges[i].get_vertices(vx);
    if (edges[i].direction == OUTGOING) {
      pr.first = vx[0]; pr.second = vx[1]; 
      edges[i].flow = edges[i].capacity - double(rgraph[pr])/10000.0;
    }
    else if (edges[i].direction == INCOMING) {
      pr.first = vx[1]; pr.second = vx[0];
      edges[i].flow = edges[i].capacity - double(rgraph[pr])/10000.0;
    }
  }
  return double(max_flow)/10000.0;
}

void Directed_Graph::compute_sinks(std::set<int>& output) const
{
  int i,j;
  bool sink;
  std::set<int> S;
  std::set<int>::const_iterator it;
  hash_map::const_iterator qt;

  output.clear();
  // A sink is a vertex all of whose edges are incoming...
  for(i=0; i<nvertex; ++i) {
    if (neighbours[i].empty()) continue;
    sink = true;
    for(it=neighbours[i].begin(); it!=neighbours[i].end(); ++it) {
      j = *it;
      S.clear();
      S.insert(i); S.insert(j);
      qt = index_table.find(S);
      if (edges[qt->second].get_direction(i,j) == OUTGOING) {
        // If this vertex has an outgoing edge, it isn't a sink...
        sink = false;
        break;
      }
    }
    if (sink) output.insert(i);
  }
}

void Directed_Graph::compute_sources(std::set<int>& output) const
{
  int i,j;
  bool source;
  std::set<int> S;
  std::set<int>::const_iterator it;
  hash_map::const_iterator qt;

  output.clear();
  // A source is a vertex all of whose edges are outgoing...
  for(i=0; i<nvertex; ++i) {
    if (neighbours[i].empty()) continue;
    source = true;
    for(it=neighbours[i].begin(); it!=neighbours[i].end(); ++it) {
      j = *it;
      S.clear();
      S.insert(i); S.insert(j);
      qt = index_table.find(S);
      if (edges[qt->second].get_direction(i,j) == INCOMING) {
        // If this vertex has an incoming edge, it isn't a source...
        source = false;
        break;
      }
    }
    if (source) output.insert(i);
  }
}
