#include "directed_graph.h"

using namespace SYNARMOSMA;

Directed_Graph::Directed_Graph() : Graph()
{

}

Directed_Graph::Directed_Graph(int n) : Graph(n)
{
  int i,j,nc = 0;
  double alpha;
  std::set<int> vx;
  Relation d;
  Random RND;

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
  Random RND;

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

int Directed_Graph::build_lattice(int n,const std::vector<int>& btype,int time)
{
  if (n < 2) throw std::invalid_argument("The number of vertices per dimension in Directed_Graph::build_lattice must be greater than one!");
  int d = (signed) btype.size();
  if (d < 2) throw std::invalid_argument("The number of dimensions in Directed_Graph::build_lattice must be greater than one!");
  if (time >= d) throw std::invalid_argument("The time dimension in Directed_Graph::build_lattice must not exceed the length of the boundary type vector!");
  for(int i=0; i<d; ++i) {
    if (btype[i] != 0 && btype[i] != 1) throw std::invalid_argument("The graph's boundary topology in in Directed_Graph::build_lattice must be either linear or toroidal!");
  }
  
  // This constructor builds a graph with the topology of a rectangular or Cartesian lattice, of 
  // dimensionality equal to the length of the vector btype and with n vertices per dimension. The 
  // btype vector determines the nature of the boundary for each dimension: 0 means it is linear, 
  // 1 means toroidal.
  int i,j,k,p,q,in1,base;
  std::set<int> S;
  std::vector<int> vx,index;
  Relation rho;

  clear();

  if (time >= 0) {
    if (btype[time] == 1) std::cout << "Warning: This directed graph will have a toroidal time dimension!" << std::endl;
  }

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
      if (j == time) {
        rho = Relation::before;
      }
      else {
        rho = Relation::disparate;
      }
      vx = index;
      if (index[j] > 0) {
        vx[j] -= 1;
        in1 = 0;
        for(k=0; k<d; ++k) {
          in1 += ipow(n,d-1-k)*vx[k];
        }
        if (i < in1) {
          add_edge(i,in1,rho);
          if (rho == Relation::before) number_directed++;
        } 
      }
      else {
        if (btype[j] == 1) {
          vx[j] = n-1;
          in1 = 0;
          for(k=0; k<d; ++k) {
            in1 += ipow(n,d-1-k)*vx[k];
          }
          if (i < in1) {
            add_edge(i,in1,rho);
            if (rho == Relation::before) number_directed++;
          } 
        }
      }
      vx = index;
      if (index[j] < n-1) {
        vx[j] += 1;
        in1 = 0;
        for(k=0; k<d; ++k) {
          in1 += ipow(n,d-1-k)*vx[k];
        }
        if (i < in1) {
          add_edge(i,in1,rho); 
          if (rho == Relation::before) number_directed++;
        }
      }
    }
  }
  return size();
}

void Directed_Graph::clear()
{
  number_directed = 0;
  nvertex = 0;
  neighbours.clear();
  edges.clear();
  index_table.clear();
}

int Directed_Graph::out_degree(int v) const
{
  if (v < 0 || v >= nvertex) throw std::invalid_argument("Illegal vertex index in the Directed_Graph::out_degree method!");
  int u,n = 0;
  std::set<int> S;
  std::set<int>::const_iterator it;
  hash_map::const_iterator qt;

  for(it=neighbours[v].begin(); it!=neighbours[v].end(); ++it) {
    u = *it;
    S.clear();
    S.insert(v); S.insert(u);
    qt = index_table.find(S);
    // Does this edge point from v to u?
    if (edges[qt->second].get_direction(v,u) == Relation::before) n++;
  }     
  return n;
}

int Directed_Graph::in_degree(int v) const
{
  if (v < 0 || v >= nvertex) throw std::invalid_argument("Illegal vertex index in the Directed_Graph::in_degree method!");
  int u,n = 0;
  std::set<int> S;
  std::set<int>::const_iterator it;
  hash_map::const_iterator qt;

  for(it=neighbours[v].begin(); it!=neighbours[v].end(); ++it) {
    u = *it;
    S.clear();
    S.insert(v); S.insert(u);
    qt = index_table.find(S);
    // Does this edge point from u to v?
    if (edges[qt->second].get_direction(v,u) == Relation::after) n++;
  }     
  return n;
}

int Directed_Graph::neutral_degree(int v) const
{
  if (v < 0 || v >= nvertex) throw std::invalid_argument("Illegal vertex index in the Directed_Graph::neutral_degree method!");
  int u,n = 0;
  std::set<int> S;
  std::set<int>::const_iterator it;
  hash_map::const_iterator qt;

  for(it=neighbours[v].begin(); it!=neighbours[v].end(); ++it) {
    u = *it;
    S.clear();
    S.insert(v); S.insert(u);
    qt = index_table.find(S);
    // Does this edge have no orientation?
    if (edges[qt->second].get_direction(v,u) == Relation::disparate) n++;
  }     
  return n;
}

void Directed_Graph::write2disk(const std::string& filename,const std::vector<std::string>& names) const
{
  int i,j;
  Relation rho;
  std::set<int> S;
  std::set<int>::const_iterator it;
  hash_map::const_iterator qt;

  std::ofstream s(filename,std::ios::trunc);

  s << "digraph G {" << std::endl;
  if (names.empty()) {
    // First all the vertices...
    for(i=0; i<nvertex; ++i) {
      s << "  \"" << 1+i << "\";" << std::endl;
    }
    // Now the edges...
    for(i=0; i<nvertex; ++i) {
      for(it=neighbours[i].begin(); it!=neighbours[i].end(); ++it) {
        j = *it;
        // We assume no multi-edges...
        if (i > j) continue;
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
        // We assume no multi-edges...
        if (i > j) continue;
        S.clear();
        S.insert(i); S.insert(j);
        qt = index_table.find(S);
        rho = edges[qt->second].get_direction(i,j);
        if (rho == Relation::before) {
          s << "  \"" << names[i] << "\" -> \"" << names[j] << "\";" << std::endl;
        }
        else if (rho == Relation::after) {
          s << "  \"" << names[j] << "\" -> \"" << names[i] << "\";" << std::endl;
        }
        else {
          s << "  \"" << names[i] << "\" -- \"" << names[j] << "\";" << std::endl;
        }
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

bool Directed_Graph::eulerian() const
{
  if (!connected()) return false;

  int i,j,n_in,n_out;
  std::set<int> S;
  std::set<int>::const_iterator jt;
  hash_map::const_iterator qt;

  for(i=0; i<nvertex; ++i) {
    n_in = 0; n_out = 0;
    for(jt=neighbours[i].begin(); jt!=neighbours[i].end(); ++jt) {
      j = *jt;
      S.clear();
      S.insert(i); S.insert(j);
      qt = index_table.find(S);
      // Does this edge point from i to j?
      if (edges[qt->second].get_direction(i,j) == Relation::before) {
        n_out++;
      }
      else if (edges[qt->second].get_direction(i,j) == Relation::after) {
        n_in++;
      } 
    }
    if (n_in != n_out) return false;
  }
  return true;
}

int Directed_Graph::compute_hamiltonian_path(bool cycle,int nattempts,std::vector<int>& path,int start) const
{
  if (!connected()) {
    path.clear();
    return -1;
  }
  int i,u,v,w,l,n,svertex,winner;
  bool visited[nvertex],solution;
  std::vector<std::pair<int,int> > candidates;
  std::set<int> S,next;
  std::set<int>::const_iterator it,jt;
  Random RND;
  hash_map::const_iterator qt;

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
        S.clear();
        S.insert(v); S.insert(u);
        qt = index_table.find(S);
        if (edges[qt->second].get_direction(v,u) != Relation::before) continue;
        n = 0;
        for(jt=neighbours[u].begin(); jt!=neighbours[u].end(); ++jt) {
          w = *jt;
          if (visited[w]) continue;
          S.clear();
          S.insert(u); S.insert(w);
          qt = index_table.find(S);
          if (edges[qt->second].get_direction(u,w) == Relation::before) n++;
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
  if (u < 0 || u >= nvertex) throw std::invalid_argument("Illegal vertex argument in the Directed_Graph::distance method!");
  if (v < 0 || v >= nvertex) throw std::invalid_argument("Illegal vertex argument in the Directed_Graph::distance method!");

  if (u == v) return 0;

  if (connected(u,v)) return 1;

  // So, we'll have to use Dijkstra's algorithm then...
  int i,n,p,md,base,l = -1;
  bool visited[nvertex];
  std::set<int> S;
  std::vector<int> distance;
  std::set<int>::const_iterator it;
  hash_map::const_iterator qt;

  for(i=0; i<nvertex; ++i) {
    distance.push_back(-1);
    visited[i] = false;
  }
  distance[u] = 0;
  do {
    base = -1;
    md = nvertex + 1;
    for(i=0; i<nvertex; ++i) {
      if (visited[i]) continue;
      n = distance[i];
      if (n == -1) continue;
      if (n < md) {
        md = n;
        base = i;
      }
    }
    if (base == -1 || base == v) break;
    visited[base] = true;
    n = distance[base] + 1;
    for(it=neighbours[base].begin(); it!=neighbours[base].end(); ++it) {
      p = *it;
      S.clear();
      S.insert(base); S.insert(p);
      qt = index_table.find(S);
      if (edges[qt->second].get_direction(base,p) != Relation::before) continue;   
      if (distance[p] < 0 || distance[p] > n) distance[p] = n;
    }
  } while(true);
  if (base == v) l = md;
  return l;
}

void Directed_Graph::compute_shortest_path(int u,int v,std::vector<int>& path) const
{
  if (u < 0 || u >= nvertex) throw std::invalid_argument("Illegal vertex argument in the Directed_Graph::compute_shortest_path method!");
  if (v < 0 || v >= nvertex) throw std::invalid_argument("Illegal vertex argument in the Directed_Graph::compute_shortest_path method!");
  if (u == v) throw std::invalid_argument("Vertex arguments are identical in the Directed_Graph::compute_shortest_path method!");

  path.clear();

  if (connected(u,v)) {
    path.push_back(v);
    return;
  }

  // So, we'll have to use Dijkstra's algorithm then...
  int i,n,p,md,base;
  bool visited[nvertex];
  std::set<int> S;
  std::vector<int> distance,antecedent;
  std::set<int>::const_iterator it;
  hash_map::const_iterator qt;

  for(i=0; i<nvertex; ++i) {
    distance.push_back(-1);
    antecedent.push_back(-1);
    visited[i] = false;
  }
  distance[u] = 0;
  do {
    base = -1;
    md = nvertex + 1;
    for(i=0; i<nvertex; ++i) {
      if (visited[i]) continue;
      n = distance[i];
      if (n == -1) continue;
      if (n < md) {
        md = n;
        base = i;
      }
    }
    if (base == -1 || base == v) break;
    visited[base] = true;
    n = distance[base] + 1;
    for(it=neighbours[base].begin(); it!=neighbours[base].end(); ++it) {
      p = *it;
      S.clear();
      S.insert(base); S.insert(p);
      qt = index_table.find(S);
      if (edges[qt->second].get_direction(base,p) != Relation::before) continue;     
      if (distance[p] < 0 || distance[p] > n) {
        distance[p] = n;
        antecedent[p] = base;
      }
    }
  } while(true);
  if (base != v) return;
  // Now calculate the reverse path...
  do {
    path.push_back(base);
    base = antecedent[base];
    if (base == u) break;
  } while(true);
  // And put it in the correct order before exiting...
  std::reverse(path.begin(),path.end());
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
  Random RND;
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
