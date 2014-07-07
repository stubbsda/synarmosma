#include "edge.h"
#include "schema.h"

#ifndef __dgraph
#define __dgraph

class Directed_Graph : public Schema {
 private:
  std::vector<Edge> edges;
  hash_map index_table;

  // Hyphantic operators
  bool foliation_x(int,int);
  bool foliation_m(int,int);
  int fission_x(int);
  int fission_m(int);
 public:
  Directed_Graph();
  Directed_Graph(int);
  virtual ~Directed_Graph();
  bool add_edge(int,int);
  bool add_edge(int,int,DIRECTION);
  int two_cycles() const;
  bool path_connected(int,int) const;
};
#endif
