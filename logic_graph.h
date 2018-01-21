#include "graph.h"
#include "propositional_system.h"

#ifndef _lgraphh
#define _lgraphh

namespace SYNARMOSMA {
  class Logic_Graph: public Graph {
   protected:  
    // The main property of this class
    Propositional_System* logic;
    // This array measures the number of atomic propositions used in a given 
    // graph neighbourhood (a vertex and its neighbours)
    std::vector<int> logical_breadth;

    double rationalize_topology();
    void compute_logical_breadth();

    virtual bool fusion(int,int);
    virtual bool foliation_x(int,int);
    virtual bool foliation_m(int,int);
    virtual int fission_x(int);
    virtual int fission_m(int);
    virtual bool add_edge(int,int);  
    virtual bool drop_vertex(int);
   public:
    Logic_Graph();
    Logic_Graph(int);
    Logic_Graph(const Logic_Graph&);
    Logic_Graph& operator =(const Logic_Graph&);
    virtual ~Logic_Graph();
    virtual void clear();
    virtual int serialize(std::ofstream&) const;
    virtual int deserialize(std::ifstream&);
    void create();
  };
}
#endif
