#include "graph.h"
#include "propositional_system.h"

#ifndef _lgraphh
#define _lgraphh

namespace SYNARMOSMA {
  /// A class representing a graph where every vertex is associated with a logical proposition.
  class Logic_Graph: public Graph {
   protected:  
    /// The main property of this class, a pointer to the Propositional_System class, 
    /// which stores the logical proposition associated with each vertex.
    Propositional_System* logic;
    /// This array measures the number of atomic propositions used in each graph 
    /// neighbourhood, i.e. a vertex and its neighbours.
    std::vector<unsigned int> logical_breadth;

    /// This method accepts a Boolean operator (by default AND) and uses propositional consistency according to it in order to determine how the graph topology should be organized; the method returns a value which measures the degree to which the graph topology diverges from this ideal. 
    double rationalize_topology(const Boolean& = Boolean::AND);
    /// This method computes the value of the vector Logic_Graph::logical_breadth and returns the total number of atomic propositions used in the entire graph.
    unsigned int compute_logical_breadth();
    /// This method first calls the Graph::fusion method on the arguments and if it returns true, then it calls the Propositional_System::drop_theorem method on the second argument and returns the value of this call.
    bool fusion(int,int) override;
    /// This method calls the Graph::fission_x method on the argument - the vertex undergoing fission - and then adds a new proposition using add_theorem(); the method returns the index of the new vertex.  
    int fission_x(int) override;
    /// This method calls the Graph::fission_m method on the argument - the vertex undergoing fission - and then adds a new proposition using add_theorem(); the method returns the index of the new vertex.  
    int fission_m(int) override;
    /// This method first calls the Graph::drop_vertex method on the argument and if it returns true, then it calls the Propositional_System::drop_theorem method on this vertex and returns the value of this call.
    bool drop_vertex(int) override;
    /// This method adds a new theorem to Logic_Graph::logic using as atomic propositions those of the argument and its neighbours.
    void add_theorem(int);
   public:
    /// The default constructor which does nothing but call the default Graph constructor.
    Logic_Graph();
    /// The principal constructor for this class, with the first argument the number of vertices in the logic graph and the second argument the percentage of these vertices that can be used as atomic propositions.
    Logic_Graph(int,double);
    /// The standard copy constructor, which copies over the properties from the source instance.
    Logic_Graph(const Logic_Graph&);
    /// The standard overloaded assignment operator, which copies over the properties from the source instance.
    Logic_Graph& operator =(const Logic_Graph&);
    /// The destructor which releases the memory from Logic_Graph::logic if Logic_Graph::nvertex is greater than zero.
    ~Logic_Graph() override;
    /// This method sets all the properties to their default value: zero vertices and empty edge and neighbour tables, as well as deleting the Logic_Graph::logic property.
    void clear() override;
    /// This method writes the instance properties to a binary disk file and returns the number of bytes written to the file.
    int serialize(std::ofstream&) const override;
    /// This method calls the clear() method on the instance and then reads the properties from a binary disk file and returns the number of bytes read.
    int deserialize(std::ifstream&) override;
  };
}
#endif
