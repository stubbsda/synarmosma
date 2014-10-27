#include "graph.h"
#include "propositional_system.h"

#ifndef _lgraphh
#define _lgraphh

namespace SYNARMOSMA {
  class Logic_Graph: public Graph {
   private:  
    Propositional_System* logic;
    // This array measures the number of atomic propositions used in a given 
    // graph neighbourhood (a vertex and its neighbours)
    std::vector<int> logical_breadth;
    int rationalize_topology();
    void compute_logical_breadth();

    virtual bool amputation(int);
    virtual bool fusion(int,int);
    virtual bool foliation_x(int,int);
    virtual bool foliation_m(int,int);
    virtual int fission_x(int);
    virtual int fission_m(int);
    virtual bool add_edge(int,int);  
   public:
    Logic_Graph();
    Logic_Graph(int);
    virtual ~Logic_Graph();
    void create();
  };
}
#endif
