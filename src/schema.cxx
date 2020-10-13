#include "schema.h"

using namespace SYNARMOSMA;

Schema::Schema()
{

}

Schema::Schema(int n)
{
  if (n <= 0) throw std::invalid_argument("The number of vertices in the Schema constructor must be greater than zero!");
  for(int i=0; i<n; ++i) {
    add_vertex();
  }
}

Schema::Schema(const Schema& source)
{
  nvertex = source.nvertex;
  neighbours = source.neighbours;
}

Schema& Schema::operator =(const Schema& source)
{
  if (this == &source) return *this;

  nvertex = source.nvertex;
  neighbours = source.neighbours;

  return *this;
}

Schema::~Schema()
{

}

void Schema::clear()
{
  nvertex = 0;
  neighbours.clear();
}

int Schema::serialize(std::ofstream& s) const
{
  int i,j,count = 0;
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
  return count;
}

int Schema::deserialize(std::ifstream& s)
{
  int i,j,k,n,count = 0;
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
  return count;
}

bool Schema::drop_edge(int n,int m)
{
  if (n < 0 || n >= nvertex) throw std::invalid_argument("A vertex argument in Schema::drop_edge does not exist!");
  if (m < 0 || m >= nvertex) throw std::invalid_argument("A vertex argument in Schema::drop_edge does not exist!");

  if (n == m) return false;

  std::set<int>::const_iterator it;

  it = neighbours[n].find(m);
  if (it == neighbours[n].end()) {
    // This edge doesn't exist...
    return false;
  }
  neighbours[n].erase(*it);
  it = neighbours[m].find(n);
  neighbours[m].erase(*it);
  return true;
}

bool Schema::add_edge(int n,int m)
{
  if (n < 0 || n >= nvertex) throw std::invalid_argument("A vertex argument in Schema::add_edge does not exist!");
  if (m < 0 || m >= nvertex) throw std::invalid_argument("A vertex argument in Schema::add_edge does not exist!");

  if (n == m) return false;

  if (neighbours[n].count(m) == 0) {
    // This edge doesn't already exist, so add it...
    neighbours[n].insert(m);
    neighbours[m].insert(n);
    return true;
  }
  return false;
}

bool Schema::consistent() const
{
  // Another utility routine, useful in debugging, that checks for various
  // pathological conditions in the neighbour table, such as vertices that have
  // a valence less than or equal to zero, or claim to be connected to themselves,
  // or have two bonds to the same vertex, or have (locally) inconsistent bonds,
  // or claim to be bonded to nonexistent vertices. If any of these conditions
  // arise, then a relatively verbose error message is printed out, and the method
  // returns false.
  int i,in1;
  std::set<int>::const_iterator it;

  if (nvertex != (signed) neighbours.size()) return false;

  // Loop through all vertices
  for(i=0; i<nvertex; ++i) {
    // Check if the valence is zero...
    if (neighbours[i].size() == 0) return false;
    // Now loop through all the neighbours
    for(it=neighbours[i].begin(); it!=neighbours[i].end(); ++it) {
      in1 = *it;
      // Check for self-bonding...
      if (in1 == i) return false;
      // Check for neighbour values that are negative or too large...
      if (in1 < 0 || in1 >= nvertex) return false;
      // Check that this vertex's neighbour list is consistent with mine... 
      if (neighbours[in1].count(i) == 0) return false;
    }
  }
  return true;
}

int Schema::distance(int v1,int v2) const
{
  if (v1 < 0 || v1 >= nvertex) throw std::invalid_argument("Illegal vertex argument in the Schema::distance method!");
  if (v2 < 0 || v2 >= nvertex) throw std::invalid_argument("Illegal vertex argument in the Schema::distance method!");

  // A method to calculate the topological distance
  // between the two vertices v1 and v2
  if (v1 == v2) return 0;

  // Handle the simplest case first, where v2 is a neighbour of v1...
  if (connected(v1,v2)) return 1;

  // So, we'll have to use Dijkstra's algorithm then...
  int i,n,md,base,l = -1;
  bool visited[nvertex];
  std::vector<int> distance;
  std::set<int>::const_iterator it;

  for(i=0; i<nvertex; ++i) {
    distance.push_back(-1);
    visited[i] = false;
  }
  distance[v1] = 0;
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
    if (base == -1 || base == v2) break;
    visited[base] = true;
    n = distance[base] + 1;
    for(it=neighbours[base].begin(); it!=neighbours[base].end(); ++it) {
      if (distance[*it] < 0 || distance[*it] > n) distance[*it] = n;
    }
  } while(true);
  if (base == v2) l = md;
  return l;
}

void Schema::compute_shortest_path(int v1,int v2,std::vector<int>& path) const
{
  if (v1 < 0 || v1 >= nvertex) throw std::invalid_argument("Illegal vertex argument in the Schema::compute_shortest_path method!");
  if (v2 < 0 || v2 >= nvertex) throw std::invalid_argument("Illegal vertex argument in the Schema::compute_shortest_path method!");
  if (v1 == v2) throw std::invalid_argument("Vertex arguments are identical in the Schema::compute_shortest_path method!");

  path.clear();

  if (connected(v1,v2)) {
    path.push_back(v2);
    return;
  }

  // So, we'll have to use Dijkstra's algorithm then...
  int i,n,v,md,base;
  bool visited[nvertex];
  std::vector<int> distance,antecedent;
  std::set<int>::const_iterator it;

  for(i=0; i<nvertex; ++i) {
    distance.push_back(-1);
    antecedent.push_back(-1);
    visited[i] = false;
  }
  distance[v1] = 0;
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
    if (base == -1 || base == v2) break;
    visited[base] = true;
    n = distance[base] + 1;
    for(it=neighbours[base].begin(); it!=neighbours[base].end(); ++it) {
      v = *it;
      if (distance[v] < 0 || distance[v] > n) {
        distance[v] = n;
        antecedent[v] = base;
      }
    }
  } while(true);
  if (base != v2) return;
  // Now calculate the reverse path...
  do {
    path.push_back(base);
    base = antecedent[base];
    if (base == v1) break;
  } while(true);
  // And put it in the correct order before exiting...
  std::reverse(path.begin(),path.end());
}

void Schema::compute_distances(pair_index& output) const
{
  int i,j,k,delta,distances[nvertex][nvertex];
  std::pair<int,int> pr;
  std::set<int>::const_iterator it;

  for(i=0; i<nvertex; ++i) {
    for(j=0; j<nvertex; ++j) {
      distances[i][j] = std::numeric_limits<int>::max();
    }
    distances[i][i] = 0;
    for(it=neighbours[i].begin(); it!=neighbours[i].end(); ++it) {
      distances[i][*it] = 1;
    }
  }

  for(k=0; k<nvertex; ++k) {
    for(i=0; i<nvertex; ++i) {
      if (distances[i][k] == std::numeric_limits<int>::max()) continue;
      for(j=0; j<nvertex; ++j) {
        if (distances[k][j] == std::numeric_limits<int>::max()) continue;
        delta = distances[i][k] + distances[k][j];
        if (delta < distances[i][j]) distances[i][j] = delta;
      }
    }
  }
  output.clear();
  for(i=0; i<nvertex; ++i) {
    for(j=1+i; j<nvertex; ++j) {
      pr.first = i; pr.second = j;
      output[pr] = distances[i][j];
    }
  }
}

bool Schema::connected() const
{
  if (nvertex == 1) return true;

  int i,n = -1,p,m = 0;
  std::vector<int> ubiquity;
  std::set<int> change,nchange;
  std::set<int>::const_iterator it,jt;

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
        nchange.insert(n);
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

int Schema::component_analysis(std::vector<int>& component) const
{
  int i,n,m,start,nc = 0;
  std::set<int> change,nchange;
  std::set<int>::const_iterator it,jt;

  component.clear();
  for(i=0; i<nvertex; ++i) {
    component.push_back(-1);
  }

  do {
    start = -1;
    for(i=0; i<nvertex; ++i) {
      if (component[i] == -1) {
        start = i;
        break;
      }
    }
    if (start == -1) break;
    change.insert(start);
    do {
      for(it=change.begin(); it!=change.end(); ++it) {
        n = *it;
        component[n] = nc;
      }
      for(it=change.begin(); it!=change.end(); ++it) {
        n = *it;
        for(jt=neighbours[n].begin(); jt!=neighbours[n].end(); ++jt) {
          m = *jt;
          if (component[m] < 0) nchange.insert(m);
        }
      }
      if (nchange.empty()) break;
      change = nchange;
      nchange.clear();
    } while(true);
    nc++;
    change.clear();
  } while(true);
  return nc;
}

int Schema::spanning_tree(std::vector<int>& tree_edges) const
{
  int ntree,n,m;
  std::set<int> current,vertices,next;
  std::set<int>::const_iterator it,jt,kt,lt;

  // A sanity check...
  if (!connected()) throw std::invalid_argument("The spanning tree calculation only makes sense for a connected schema!");

  current.insert(0);

  do {
    for(it=current.begin(); it!=current.end(); ++it) {
      vertices.insert(*it);
    }
    for(it=current.begin(); it!=current.end(); ++it) {
      n = *it;
      for(jt=neighbours[n].begin(); jt!=neighbours[n].end(); ++jt) {
        m = *jt;
        kt = std::find(vertices.begin(),vertices.end(),m);
        if (kt != vertices.end()) continue;
        lt = std::find(next.begin(),next.end(),m);
        if (lt != next.end()) continue;
        next.insert(m);
        // Now add this edge to the spanning tree...
        if (n < m) {
          tree_edges.push_back(n); tree_edges.push_back(m);
        }
        else {
          tree_edges.push_back(m); tree_edges.push_back(n);
        }
      }
    }
    if (next.empty()) break;
    current = next;
    next.clear();
  } while(true);
  ntree = (signed) tree_edges.size();

  if (nvertex != (1 + ntree/2)) throw std::runtime_error("The size of the spanning tree is inconsistent with the schema's size!");

  return ntree;
}



