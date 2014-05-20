#include "directed_graph.h"

extern Random RND;

Directed_Graph::Directed_Graph() : Schema() 
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
  name += (d == FORWARD) ? "1" : "-1";
  qt = index_table.find(name);
  if (qt != index_table.end()) return false;
  edges.push_back(Edge(v1,v2,d));
  index_table[name] = (signed) edges.size() - 1;
  return true;
}

