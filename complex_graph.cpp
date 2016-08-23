#include "complex_graph.h"

using namespace SYNARMOSMA;

Complex_Graph::Complex_Graph()
{
  nvertex = 0;
}

Complex_Graph::Complex_Graph(int order)
{
  assert(order > 0);
  nvertex = order;
  neighbours = new std::vector<int>[nvertex];
}

Complex_Graph::Complex_Graph(const Complex_Graph& source)
{
  if (nvertex > 0) delete[] neighbours;
  nvertex = source.nvertex;
  neighbours = new std::vector<int>[nvertex];
  for(int i=0; i<nvertex; ++i) {
    neighbours[i] = source.neighbours[i];
  }
}

Complex_Graph& Complex_Graph::operator =(const Complex_Graph& source)
{
  if (this == &source) return *this;
  if (nvertex > 0) delete[] neighbours;
  nvertex = source.nvertex;
  neighbours = new std::vector<int>[nvertex];
  for(int i=0; i<nvertex; ++i) {
    neighbours[i] = source.neighbours[i];
  }
  return *this;
}
 
Complex_Graph::~Complex_Graph()
{
  if (nvertex > 0) delete[] neighbours;
}

void Complex_Graph::add_edge(int u,int v)
{
  assert(u >= 0 && u < nvertex);
  assert(v >= 0 && v < nvertex);
  neighbours[u].push_back(v);
  neighbours[v].push_back(u);
}

void Complex_Graph::display() const
{
  int i,j,n,k = 0;
  for(i=0; i<nvertex; ++i) {
    n = (signed) neighbours[i].size();
    for(j=0; j<n; ++j) {
      std::cout << i << ":" << neighbours[i][j] << std::endl;
      k++;
    }
  }
  std::cout << "This graph has " << nvertex << " vertices and " << k/2 << " edges." << std::endl;
}

int Complex_Graph::DFS_bridge(int u,int v,int dcount,int* low,int* pre,hash_map& bridge_index) const
{
  int w,output = 0,dc = dcount + 1;
  std::vector<int>::const_iterator it;

  pre[v] = dc;
  low[v] = pre[v];
  for(it=neighbours[v].begin(); it!=neighbours[v].end(); ++it) {
    w = *it;
    if (w == v) continue;
    if (pre[w] == -1) {
      output += DFS_bridge(v,w,dc,low,pre,bridge_index);
      low[v] = std::min(low[v],low[w]); 
      if (low[w] == pre[w]) {
        //std::cout << "Found bridge for edge " << v << ":" << w << std::endl;
        std::set<int> S;
        S.insert(v); S.insert(w);
        bridge_index[S] = output;
        output++;
      }
    }
    else if (w != u) {
      low[v] = std::min(low[v],pre[w]);
    }
  }
  return output;
}

int Complex_Graph::compute_bridges(hash_map& bridge_index) const
{
  int i,bcount = 0;
  int low[nvertex],pre[nvertex];

  for(i=0; i<nvertex; ++i) {
    low[i] = -1;
    pre[i] = -1;
  }
  bridge_index.clear();
  for(i=0; i<nvertex; ++i) {
    if (pre[i] == -1) bcount += DFS_bridge(i,i,0,low,pre,bridge_index);
  }
  assert(bcount == bridge_index.size());
  return bcount;
}

int Complex_Graph::get_candidates(std::vector<int>& output) const
{
  // The return value is the number of bridges in this graph, while the 
  // "output" vector is the list of potential edges for contraction or 
  // deletion, listed by pairs of vertices.
  int i,j,k,l,n,nb;
  bool first;
  std::set<int> S;
  hash_map bridge_index;
  hash_map::const_iterator qt;

  nb = compute_bridges(bridge_index);
  output.clear();
  for(i=0; i<nvertex; ++i) {
    n = (signed) neighbours[i].size();
    first = true;
    for(j=0; j<n; ++j) {
      k = neighbours[i][j];
      if (k <= i) continue;
      // Check if this is a bridge
      S.insert(i); S.insert(k);
      qt = bridge_index.find(S);
      if (qt == bridge_index.end()) {
        output.push_back(i); output.push_back(k);
      }
      else {
        // Since this is a multi-graph, we need to check for overcounting of bridges...
        l = std::count(neighbours[i].begin(),neighbours[i].end(),k);
        if (l > 1) {
          if (first) {
            nb -= 1;
            first = false;
          }
          output.push_back(i); output.push_back(k);
        }
      }
      S.clear();
    }
  }
  return nb;
}

void Complex_Graph::contract(int u,int v,Complex_Graph* output) const
{
  // This involves fusing together two vertices, so the output graph will have 
  // one less edge and one less vertex, which will mean re-indexing all of the 
  // neighbour vectors.
  int i,j,n,vmin,vmax; 
  bool first = true;
  std::vector<int> nvector;

  if (u > v) {
    vmin = v;
    vmax = u;
  }
  else {
    vmin = u;
    vmax = v;
  }
  // We will fuse vmax into vmin
  for(i=0; i<nvertex; ++i) {
    if (i == vmax) continue;
    n = (signed) neighbours[i].size();
    for(j=0; j<n; ++j) {
      if (neighbours[i][j] == vmax) {
        if (i == vmin) {
          if (first) {
            first = false;
          }
          else {
            nvector.push_back(vmin);
          }
        }
        else {
          nvector.push_back(vmin);
        }
      }
      else if (neighbours[i][j] > vmax) {
        nvector.push_back(neighbours[i][j] - 1);
      }
      else {
        nvector.push_back(neighbours[i][j]);
      }
    }
    output->neighbours[i] = nvector;
    nvector.clear();
  }
  n = (signed) neighbours[vmax].size();
  for(i=0; i<n; ++i) {
    j = neighbours[vmax][i];
    if (j == vmin) {
      continue;
    }
    else if (j == vmax) {
      output->neighbours[vmin].push_back(vmin);
    }
    else if (j > vmax) {
      output->neighbours[vmin].push_back(j-1);
    }
    else {
      output->neighbours[vmin].push_back(j);
    }
  }
  for(i=vmax; i<nvertex-1; ++i) {
    output->neighbours[i] = output->neighbours[i+1]; 
  }
  output->nvertex -= 1;
}

void Complex_Graph::remove(int u,int v,Complex_Graph* output) const
{
  // The much simpler case of removing an edge...
  int i,n = (signed) output->neighbours[u].size();
  bool first = true;
  std::vector<int> t;
  for(i=0; i<n; ++i) {
    if (first && output->neighbours[u][i] == v) {
      first = false;
      continue;
    }
    t.push_back(output->neighbours[u][i]);
  }
  output->neighbours[u] = t;
  t.clear();
  first = true;
  n = (signed) output->neighbours[v].size();
  for(i=0; i<n; ++i) {
    if (first && output->neighbours[v][i] == u) {
      first = false;
      continue;
    }
    t.push_back(output->neighbours[v][i]);
  }
  output->neighbours[v] = t;  
}
