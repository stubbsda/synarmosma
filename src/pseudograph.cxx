#include "pseudograph.h"

using namespace SYNARMOSMA;

Pseudograph::Pseudograph()
{

}

Pseudograph::Pseudograph(int order)
{
  if (order <= 0) throw std::invalid_argument("The order of a pseudograph must be greater than zero!");

  std::set<int> empty;
  nvertex = order;
  for(int i=0; i<nvertex; ++i) {
    entourage.push_back(std::pair<int,std::set<int> >(0,empty));
  }
}

Pseudograph::Pseudograph(const Pseudograph& source)
{
  nvertex = source.nvertex;
  entourage = source.entourage;
  edges = source.edges;
}

Pseudograph& Pseudograph::operator =(const Pseudograph& source)
{
  if (this == &source) return *this;

  nvertex = source.nvertex;
  entourage = source.entourage;
  edges = source.edges;

  return *this;
}

Pseudograph::~Pseudograph()
{

}

void Pseudograph::clear()
{
  nvertex = 0;
  entourage.clear();
  edges.clear();
}

int Pseudograph::serialize(std::ofstream& s) const
{
  int i,n,p,count = 0;
  std::set<int> T;
  std::set<int>::const_iterator it;

  s.write((char*)(&nvertex),sizeof(int)); count += sizeof(int);
  for(i=0; i<nvertex; ++i) {
    n = entourage[i].first;
    s.write((char*)(&n),sizeof(int)); count += sizeof(int);
    T = entourage[i].second;
    n = (signed) T.size();
    s.write((char*)(&n),sizeof(int)); count += sizeof(int);
    for(it=T.begin(); it!=T.end(); ++it) {
      p = *it;
      s.write((char*)(&p),sizeof(int)); count += sizeof(int);
    }
  }

  n = (signed) edges.size();
  for(i=0; i<n; ++i) {
    count += edges[i].serialize(s);
  }

  return count;
}

int Pseudograph::deserialize(std::ifstream& s)
{
  int i,j,n,p,count = 0;
  std::set<int> T;
  std::pair<int,std::set<int> > pr;
  Edge e;

  clear();

  s.read((char*)(&nvertex),sizeof(int)); count += sizeof(int);
  for(i=0; i<nvertex; ++i) {
    s.read((char*)(&n),sizeof(int)); count += sizeof(int);
    pr.first = n;
    s.read((char*)(&n),sizeof(int)); count += sizeof(int);
    T.clear();
    for(j=0; j<n; ++j) {
      s.read((char*)(&p),sizeof(int)); count += sizeof(int);
      T.insert(p);
    }
    pr.second = T;
    entourage.push_back(pr);
  }

  s.read((char*)(&n),sizeof(int)); count += sizeof(int);
  for(i=0; i<n; ++i) {
    count += e.deserialize(s);
    edges.push_back(e);
  }

  if (!consistent()) throw std::runtime_error("The Pseudograph instance read from disk is inconsistent!");

  return count;
}

bool Pseudograph::consistent() const
{
  if (nvertex < 0) return false;
  int i,n,vx[2];
  std::set<int> S;
  std::set<int>::const_iterator it;
  const int ne = (signed) edges.size();

  for(i=0; i<nvertex; ++i) {
    if (entourage[i].first < 0) return false;
    S = entourage[i].second;
    for(it=S.begin(); it!=S.end(); ++it) {
      n = *it;
      if (n < 0) return false;
      if (n >= ne) return false;
      edges[n].get_vertices(vx);
      if (vx[0] < 0 || vx[0] >= nvertex) return false;
      if (vx[1] < 0 || vx[1] >= nvertex) return false;
      if (vx[0] != i && vx[1] != i) return false;
    }
  }
  for(i=0; i<ne; ++i) {
    edges[i].get_vertices(vx);
    if (vx[0] < 0 || vx[0] >= nvertex) return false;
    if (vx[1] < 0 || vx[1] >= nvertex) return false;
    if (entourage[vx[0]].second.count(i) == 0) return false;
    if (entourage[vx[1]].second.count(i) == 0) return false;
  }

  return true;
}

int Pseudograph::in_degree(int v) const
{
  int i,u,vx[2],n = 0;
  std::set<int> S = entourage[v].second;
  std::set<int>::const_iterator it;

  for(it=S.begin(); it!=S.end(); ++it) {
    i = *it;
    if (edges[i].get_direction() == Relation::disparate) continue;
    edges[i].get_vertices(vx);
    u = (vx[0] == v) ? vx[1] : vx[0];
    if (edges[i].get_direction(v,u) == Relation::after) n++;
  }
  return n;
}

int Pseudograph::out_degree(int v) const
{
  int i,u,vx[2],n = 0;
  std::set<int> S = entourage[v].second;
  std::set<int>::const_iterator it;

  for(it=S.begin(); it!=S.end(); ++it) {
    i = *it;
    if (edges[i].get_direction() == Relation::disparate) continue;
    edges[i].get_vertices(vx);
    u = (vx[0] == v) ? vx[1] : vx[0];
    if (edges[i].get_direction(v,u) == Relation::before) n++;
  }
  return n;
}

int Pseudograph::multi_degree(int u,int v) const
{
  int w,vx[2],n = 0;
  std::set<int> S = entourage[u].second;
  std::set<int>::const_iterator it;

  for(it=S.begin(); it!=S.end(); ++it) {
    edges[*it].get_vertices(vx);
    w = (vx[0] == u) ? vx[1] : vx[0];
    if (w == v) n++;
  }

  return n;
}

int Pseudograph::multi_degree(int u,int v,Relation d) const
{
  int w,vx[2],n = 0;
  std::set<int> S = entourage[u].second;
  std::set<int>::const_iterator it;

  for(it=S.begin(); it!=S.end(); ++it) {
    edges[*it].get_vertices(vx);
    w = (vx[0] == u) ? vx[1] : vx[0];
    if (w == v) {
      if (edges[*it].get_direction() == d) n++;
    }
  }

  return n;
}

void Pseudograph::write2disk(const std::string& filename,const std::vector<std::string>& names) const
{
  int i,j,n,vx[2];
  std::set<int> S;
  std::set<int>::const_iterator it;
  Relation d;

  std::ofstream s(filename,std::ios::trunc);

  s << "digraph G {" << std::endl;
  if (names.empty()) {
    // First all the vertices...
    for(i=0; i<nvertex; ++i) {
      s << "  \"" << 1+i << "\";" << std::endl;
    }
    // Now the edges...
    for(i=0; i<nvertex; ++i) {
      n = entourage[i].first;
      for(j=0; j<n; ++j) {
        s << "  \"" << 1+i << "\" -- \"" << 1+i << "\";" << std::endl;
      }
      S = entourage[i].second;
      for(it=S.begin(); it!=S.end(); ++it) {
        edges[*it].get_vertices(vx);
        if (vx[0] > vx[1]) continue;
        d = edges[*it].get_direction();
        if (d == Relation::disparate) {
          s << "  \"" << 1+vx[0] << "\" -- \"" << 1+vx[1] << "\";" << std::endl;
        }
        else {
          if (d == Relation::before) {
            s << "  \"" << 1+vx[0] << "\" -> \"" << 1+vx[1] << "\";" << std::endl;
          }
          else {
            s << "  \"" << 1+vx[1] << "\" -> \"" << 1+vx[0] << "\";" << std::endl;
          }
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
      n = entourage[i].first;
      for(j=0; j<n; ++j) {
        s << "  \"" << names[i] << "\" -- \"" << names[i] << "\";" << std::endl;
      }
      S = entourage[i].second;
      for(it=S.begin(); it!=S.end(); ++it) {
        edges[*it].get_vertices(vx);
        if (vx[0] > vx[1]) continue;
        d = edges[*it].get_direction();
        if (d == Relation::disparate) {
          s << "  \"" << names[vx[0]] << "\" -- \"" << names[vx[1]] << "\";" << std::endl;
        }
        else {
          if (d == Relation::before) {
            s << "  \"" << names[vx[0]] << "\" -> \"" << names[vx[1]] << "\";" << std::endl;
          }
          else {
            s << "  \"" << names[vx[1]] << "\" -> \"" << names[vx[0]] << "\";" << std::endl;
          }
        }
      }
    }
  }
  s << "}" << std::endl;
  s.close();
}

void Pseudograph::add_edge(int u,int v,Relation d)
{
  if (u < 0 || u >= nvertex) throw std::invalid_argument("Illegal vertex index in the Pseudograph::add_edge method!");
  if (v < 0 || v >= nvertex) throw std::invalid_argument("Illegal vertex index in the Pseudograph::add_edge method!");
  // The simple case of a self-loop...
  if (u == v) {
    entourage[u].first += 1;
    return;
  }
  // The more complex case of an edge between two distinct vertices...
  int n = (signed) edges.size();
  edges.push_back(Edge(u,v,0.0,d));
  entourage[u].second.insert(n);
  entourage[v].second.insert(n);
}

int Pseudograph::DFS_bridge(int u,int v,int dcount,int* low,int* pre,hash_map& bridge_index) const
{
  int w,vx[2],output = 0,dc = dcount + 1;
  std::set<int> S = entourage[v].second;
  std::set<int>::const_iterator it;

  pre[v] = dc;
  low[v] = pre[v];
  for(it=S.begin(); it!=S.end(); ++it) {
    // We ignore directionality here...
    edges[*it].get_vertices(vx);
    w = (vx[0] == v) ? vx[1] : vx[0];
    if (pre[w] == -1) {
      output += DFS_bridge(v,w,dc,low,pre,bridge_index);
      low[v] = std::min(low[v],low[w]);
      if (low[w] == pre[w]) {
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

int Pseudograph::compute_bridges(hash_map& bridge_index) const
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
  return bcount;
}

int Pseudograph::get_candidates(std::vector<int>& output) const
{
  // The return value is the number of bridges in this graph, while the
  // "output" vector is the list of potential edges for contraction or
  // deletion, listed by pairs of vertices.
  int i,j,nb,vx[2];
  bool first;
  std::set<int> S,T;
  std::set<int>::const_iterator it;
  hash_map bridge_index;
  hash_map::const_iterator qt;

  nb = compute_bridges(bridge_index);
  output.clear();
  for(i=0; i<nvertex; ++i) {
    first = true;
    S = entourage[i].second;
    for(it=S.begin(); it!=S.end(); ++it) {
      edges[*it].get_vertices(vx);
      j = (vx[0] == i) ? vx[1] : vx[0];
      if (j < i) continue;
      // Check if this is a bridge
      T.insert(i); T.insert(j);
      qt = bridge_index.find(T);
      if (qt == bridge_index.end()) {
        output.push_back(i); output.push_back(j);
      }
      else {
        // Since this is a multi-graph, we need to check for overcounting of bridges...
        if (multi_degree(i,j) > 1) {
          if (first) {
            nb -= 1;
            first = false;
          }
          output.push_back(i); output.push_back(j);
        }
      }
      T.clear();
    }
  }
  return nb;
}

bool Pseudograph::contract(int u,int v,Pseudograph* output) const
{
  // This involves fusing together two vertices, so the output graph will have
  // one less edge and one less vertex, which will mean re-indexing all of the
  // neighbour vectors.
  if (u < 0 || u >= nvertex) throw std::invalid_argument("Illegal vertex index in the Pseudograph::contract method!");
  if (v < 0 || v >= nvertex) throw std::invalid_argument("Illegal vertex index in the Pseudograph::contract method!");
  if (u == v) throw std::invalid_argument("The vertex indices must be distinct in the Pseudograph::contract method!");
  int i,c,vmin,vmax,vx[2],ne,nv = output->nvertex;
  std::set<int> S,eliminate;
  std::set<int>::const_iterator it;
  std::set<int>::reverse_iterator rit;

  if (u > v) {
    vmin = v;
    vmax = u;
  }
  else {
    vmin = u;
    vmax = v;
  }
  // We will fuse vmax into vmin
  S = output->entourage[vmax].second;
  for(it=S.begin(); it!=S.end(); ++it) {
    i = *it;
    output->edges[i].get_vertices(vx);
    if (vx[0] == vmin || vx[1] == vmin) eliminate.insert(i);
  }
  if (eliminate.empty()) return false;

  // Add the self-loops of vmax to vmin and all but one of the vmin:vmax edges become 
  // self-loops of vmin...
  output->entourage[vmin].first += (eliminate.size() - 1) + output->entourage[vmax].first;

  // Now remove the edge(s)...
  for(rit=eliminate.rbegin(); rit!=eliminate.rend(); ++rit) {
    c = *rit;
    output->edges.erase(output->edges.begin() + c);
  }

  // Next re-index all of the vertex indices in the edges and remove the vertex...
  ne = (signed) output->edges.size();
  for(i=0; i<ne; ++i) {
    output->edges[i].get_vertices(vx);
    if (vx[0] > vmax) {
      vx[0] -= 1;
    }
    else if (vx[0] == vmax) {
      vx[0] = vmin;
    }
    if (vx[1] > vmax) {
      vx[1] -= 1;
    }
    else if (vx[1] == vmax) {
      vx[1] = vmin;
    }
    output->edges[i].set_vertices(vx[0],vx[1]);
  }
  for(i=vmax; i<nv-1; ++i) {
    output->entourage[i] = output->entourage[i+1];
  }
  output->nvertex -= 1;

  // Finally fix the vertex entourages...
  for(i=0; i<output->nvertex; ++i) {
    output->entourage[i].second.clear();
  }
  for(i=0; i<ne; ++i) {
    output->edges[i].get_vertices(vx);
    output->entourage[vx[0]].second.insert(i);
    output->entourage[vx[1]].second.insert(i);
  }

  return true;
}

bool Pseudograph::remove(int u,int v,Pseudograph* output) const
{
  if (u < 0 || u >= nvertex) throw std::invalid_argument("Illegal value for the vertex index in Pseudograph::remove!");
  if (v < 0 || v >= nvertex) throw std::invalid_argument("Illegal value for the vertex index in Pseudograph::remove!");

  int i,j,candidate = -1,vx[2],nv = output->nvertex;
  std::set<int> T,S = output->entourage[u].second;
  std::set<int>::const_iterator it;

  for(it=S.begin(); it!=S.end(); ++it) {
    i = *it;
    if (edges[i].get_direction() != Relation::disparate) continue;
    edges[i].get_vertices(vx);
    if (vx[0] == v || vx[1] == v) {
      candidate = i;
      break;
    }
  }
  // Impossible to find an undirected edge connecting these two vertices...
  if (candidate == -1) return false;

  S.erase(candidate);
  output->entourage[u].second = S;

  S = output->entourage[v].second;
  S.erase(candidate);
  output->entourage[v].second = S;

  // Now remove the edge and re-index all of the entourages...
  for(i=0; i<nv; ++i) {
    T.clear();
    S = output->entourage[i].second;
    for(it=S.begin(); it!=S.end(); ++it) {
      j = *it;
      if (j > candidate) {
        T.insert(j - 1);
      }
      else {
        T.insert(j);
      }
    }
    output->entourage[i].second = T;
  }
  output->edges.erase(output->edges.begin() + candidate);

  return true;
}

