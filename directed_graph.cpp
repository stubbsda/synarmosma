#include "directed_graph.h"

extern Random RND;

Directed_Graph::Directed_Graph() : Schema() 
{

}

Directed_Graph::Directed_Graph(int n) : Schema(n)
{

}

Directed_Graph::~Directed_Graph()
{

}

bool Directed_Graph::add_edge(int v1,int v2)
{
  bool frwd,bkwd;
  std::string name,base = make_key(v1,v2);
  hash_map::const_iterator qt;

  name = base + ":1";
  qt = index_table.find(name);
  frwd = (qt != index_table.end()) ? true : false;

  name = base + ":-1";
  qt = index_table.find(name);
  bkwd = (qt != index_table.end()) ? true : false;  

  if (frwd && bkwd) return false;

  neighbours[v1].insert(v2);
  neighbours[v2].insert(v1);

  if (frwd) {
    // Add the backward-oriented edge
    edges.push_back(Edge(v1,v2,BACKWARD));
    index_table[base+":-1"] = (signed) edges.size() - 1;
  }
  else if (bkwd) {
    // Add the forward-oriented edge
    edges.push_back(Edge(v1,v2,FORWARD));
    index_table[base+":1"] = (signed) edges.size() - 1;
  }
  else {
    // Randomly choose a direction...
    if (RND.irandom(2) == 0) {
      edges.push_back(Edge(v1,v2,FORWARD));
      index_table[base+":1"] = (signed) edges.size() - 1;
    }
    else {
      edges.push_back(Edge(v1,v2,BACKWARD));
      index_table[base+":-1"] = (signed) edges.size() - 1;
    }
  }
  return true;
}

bool Directed_Graph::add_edge(int v1,int v2,DIRECTION d)
{
  hash_map::const_iterator qt;
  std::string name = make_key(v1,v2) + ":";
  if (v1 < v2) {
    name += (d == FORWARD) ? "1" : "-1";
  }
  else {
    name += (d == FORWARD) ? "-1" : "1";
  }
  qt = index_table.find(name);
  if (qt != index_table.end()) return false;
  neighbours[v1].insert(v2);
  neighbours[v2].insert(v1);
  edges.push_back(Edge(v1,v2,d));
  index_table[name] = (signed) edges.size() - 1;
  return true;
}

bool Directed_Graph::path_connected(int u,int v) const
{
  // Is it possible to get from u to v following the orientation of the graph edges?
  int i,j;
  bool output = false;
  std::set<int> current,next;
  std::set<int>::const_iterator it,jt;
  hash_map::const_iterator qt;

  current.insert(u);
  do {
    for(it=current.begin(); it!=current.end(); ++it) {
      i = *it;
      for(jt=neighbours[i].begin(); jt!=neighbours[i].end(); ++jt) {
        j = *jt;
        if (j < i) continue;
        qt = index_table.find(make_key(i,j) + ":1");
        if (qt == index_table.end()) continue;
        next.insert(j);
      }
    }
    if (next.empty()) break;
    if (next.count(v) > 0) {
      output = true;
      break;
    }
    current = next;
    next.clear();
  } while(true);
  return output;
}

int Directed_Graph::two_cycles() const
{
  // Return how many two-cycles there are in this directed graph
  int i,j,t_cycle = 0;
  std::string key;
  std::set<int>::const_iterator it;
  hash_map::const_iterator qt;

  for(i=0; i<nvertex; ++i) {
    for(it=neighbours[i].begin(); it!=neighbours[i].end(); ++it) {
      j = *it;
      if (j < i) continue;
      key = make_key(i,j);
      qt = index_table.find(key + ":1");
      if (qt == index_table.end()) continue;
      // Now see if there is an edge going back from j to i...
      qt = index_table.find(key + ":-1");
      if (qt == index_table.end()) continue;
      t_cycle++;
    }
  }
  return t_cycle;
}

