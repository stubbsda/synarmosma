#include "schema.h"

using namespace SYNARMOSMA;

Schema::Schema()
{
  nvertex = 0;
}

Schema::Schema(int n)
{
  nvertex = 0;
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

void Schema::serialize(std::ofstream& s) const
{
  int i,j;
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
}

void Schema::deserialize(std::ifstream& s)
{
  int i,j,k,n;
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
}

bool Schema::drop_edge(int n,int m)
{
#ifdef DEBUG
  assert(n >= 0 && n < nvertex);
  assert(m >= 0 && m < nvertex);
#endif
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
#ifdef DEBUG
  assert(n >= 0 && n < nvertex);
  assert(m >= 0 && m < nvertex);
#endif
  if (n == m) return false;

  if (neighbours[n].count(m) == 0) {
    // This edge doesn't already exist, so add it...
    neighbours[n].insert(m);
    neighbours[m].insert(n);
    return true;
  }
  return false;
}

bool Schema::positive_valence() const
{
  // This method just checks if there are any vertices with no connections,
  // in which case it returns false, true otherwise
  for(int i=0; i<nvertex; ++i) {
    if (neighbours[i].empty()) return false;
  }
  return true;
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

  // Loop through all vertices
  for(i=0; i<nvertex; ++i) {
    // Check if the valence is weird
    if (neighbours[i].size() == 0) return false;
    // Now loop through all the neighbours
    for(it=neighbours[i].begin(); it!=neighbours[i].end(); ++it) {
      in1 = *it;
      // Check for self-bonding
      if (in1 == i) return false;
      // Check for neighbour values that are negative or too large
      if (in1 >= nvertex) return false;
      // If this neighbour vertex is local, then check to see if it too lists
      // i as a neighbour
      if (neighbours[in1].count(i) != 1) {
        std::cerr << "Schema inconsistent at " << i << " and " << in1 << std::endl;
        std::exit(1);
      }
    }
  }
  return true;
}

void Schema::components(std::vector<int>& csize,std::vector<int>& celements) const
{
  int i,start,nc,noe;
  std::vector<int> component;
  std::set<int> change,nchange;
  std::set<int>::const_iterator it1,it2;

  nc = 0;
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
    noe = 0;
    do {
      for(it1=change.begin(); it1!=change.end(); it1++) {
        component[*it1] = nc;
        celements.push_back(*it1);
      }
      noe += change.size();
      for(it1=change.begin(); it1!=change.end(); it1++) {
        i = *it1;
        for(it2=neighbours[i].begin(); it2!=neighbours[i].end(); it2++) {
          if (component[*it2] < 0) nchange.insert(*it2);
        }
      }
      if (nchange.empty()) break;
      change = nchange;
      nchange.clear();
    } while(true);
    csize.push_back(noe);
    nc++;
    change.clear();
  } while(true);
}

int Schema::distance(int v1,int v2) const
{
  // A method to calculate the topological distance
  // between the two vertices v1 and v2
  if (v1 == v2) return 0;

  // Handle the simplest case first, where v2 is a neighbour of v1...
  if (neighbours[v1].count(v2) == 1) return 1;

  // So, we'll have to use Dijkstra's algorithm then...
  int i,n,v,md,base,l = -1;
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
    if (base == -1) break;
    if (base == v2) {
      l = md;
      break;
    }
    visited[base] = true;
    for(it=neighbours[base].begin(); it!=neighbours[base].end(); ++it) {
      v = *it;
      n = distance[base] + 1;
      if (distance[v] < 0 || distance[v] > n) distance[v] = n;
    }
  } while(true);
  return l;
}

void Schema::compute_distances(edge_hash& output) const
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
#ifdef DEBUG
  assert(connected());
#endif

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
#ifdef DEBUG
  assert(nvertex == (1+ntree/2));
#endif
  return ntree;
}



