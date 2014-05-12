#include "graph.h"
#include "propositional_system.h"

#ifndef _lgraphh
#define _lgraphh

class Logic_Graph: public Graph {
 private:  
  Propositional_System* logic;
  // This array measures the number of atomic propositions used in a given 
  // graph neighbourhood (a vertex and its neighbours)
  int* logical_breadth;
  void rationalize_topology();
  void compute_logical_breadth();
  
 public:
  Logic_Graph();
  Logic_Graph(int);
  virtual ~Logic_Graph();
  void create();
};
#endif
