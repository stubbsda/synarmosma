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
    bool fusion(int,int) override;
    int fission_x(int) override;
    int fission_m(int) override;
    bool drop_vertex(int) override;
   public:
    Logic_Graph();
    Logic_Graph(int);
    Logic_Graph(const Logic_Graph&);
    Logic_Graph& operator =(const Logic_Graph&);
    ~Logic_Graph() override;
    void clear() override;
    int serialize(std::ofstream&) const override;
    int deserialize(std::ifstream&) override;
    void create();
  };
}
#endif
