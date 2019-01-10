#include "graph.h"

#ifndef __dgraph
#define __dgraph

namespace SYNARMOSMA {
  /// A class representing a (partially) directed graph, i.e. a graph in which some or all of the edges have an orientation.
  class Directed_Graph: public Graph {
   protected:
    /// This non-negative integer represents the number of graph edges which have a 
    /// direction, i.e. equal to Relation::before or Relation::after. 
    unsigned int number_directed = 0;
    /// This recursive method is used to compute whether or not this directed graph contains a cycle by constructing all possible outgoing paths of successively longer length and testing if the path has come back to its starting vertex, returning true if this is so.
    bool directed_cycle(const std::vector<int>&,int,int) const;

   public:
    /// The default constructor which does nothing.
    Directed_Graph();
    /// This constructor builds a complete directed graph with the number of vertices specified by the argument and the edges are undirected with 70% probability, forward pointing with 15% probability and similarly for backward pointing. 
    Directed_Graph(int);
    /// This constructor builds a directed graph with the number of vertices specified by the argument and the probability of an edge connecting a pair of vertices determined by the second argument. The edge orientation has a 70% chance of being undirected, 15% forward pointing and 15% backward pointing.
    Directed_Graph(int,double);
    /// The standard copy constructor which copies over the inherited Graph properties as well as Directed_Graph::number_directed from the source instance.
    Directed_Graph(const Directed_Graph&);
    /// The standard overloaded assignment which copies over the inherited Graph properties as well as Directed_Graph::number_directed from the source instance.
    Directed_Graph& operator =(const Directed_Graph&);
    /// The destructor, which for this class does nothing. 
    ~Directed_Graph() override;
    /// This method writes the edge's properties to a binary disk file and returns the number of bytes written to the file.
    int serialize(std::ofstream&) const override;
    /// This method calls the clear() method on the instance and then reads the properties from a binary disk file and returns the number of bytes read.
    int deserialize(std::ifstream&) override;
    /// This method adds an edge between the vertices specified by the first two arguments and exists to ensure that the appropriate method is called for this class; the method returns false if the edge already exists. 
    virtual inline bool add_edge(int u,int v) override {return add_edge(u,v,Relation::disparate,0.0);};
    /// This method adds an edge between the vertices specified by the first two arguments and the edge length; it exists to ensure that the appropriate method is called for this class and returns false if the edge already exists.
    virtual inline bool add_edge(int u,int v,double kappa) override {return add_edge(u,v,Relation::disparate,kappa);};
    /// This method adds an edge between the vertices specified by the first two arguments, followed optionally by the orientation and length of the edge. The method returns false if the edge already exists. 
    virtual bool add_edge(int,int,Relation,double = 0.0);
    /// This method to compute the maximum network flow from a source vertex to a sink vertex, respecting the orientation of the graph's edges.
    double compute_flow(int,int) override;
    /// This method computes the oriented distance between the two vertices given as arguments, i.e. the length of the shortest path connecting the two vertices that only uses directed edges, returning -1 if no such path exists.
    int distance(int,int) const override;
    /// This method computes the complete set of combinatorial distances (respecting the edge orientation) in the directed graph and stores the result as an unordered map linking pairs of vertices (i,j) where i < j and the distance d between them. 
    void compute_distances(pair_index&) const override; 
    /// This method returns true if every vertex can reach every other vertex using an oriented path, false otherwise.
    bool connected() const override;
    /// This method changes the orientation of the edge connecting the two vertices that are the method's arguments; if the edge exists, the method returns true and modifies the direction according to the following algorithm. If the edge is undirected, a coin toss determines where it is forward or backward pointing; if the edge is forward or backward pointing, with 25% probability it is reversed and with 75% probability made undirected.
    bool mutate_edge(int,int);
    /// This method computes the number of directed edges in the graph, i.e. the value of Directed_Graph::number_directed.
    inline void compute_directedness();
    /// This method returns the value of the Directed_Graph::number_directed property.
    inline unsigned int directedness() const {return number_directed;};
    /// This method returns true when this instance is a DAG (directed acyclic graph), i.e. every edge is directed and there are no cycles.
    bool DAG() const {return (acyclic() && number_directed == (unsigned) size());};
    /// This method returns true if this graph contains no directed cycles, i.e. a set of edges which allow to return to the point of depart while respecting the edge orientation. 
    bool acyclic() const;
    /// This method returns true if between every pair of vertices \f$(x,y)\f$ in this graph, there is at most one directed path leading from \f$x\f$ to \f$y\f$.  
    bool singly_connected() const;
    /// This method computes for each vertex \f$x\f$ the number of incoming edges \f$y\to x\f$ and returns the maximum of this value across all vertices in the directed graph. 
    int maximum_parents() const;
    /// This method returns true if there is a directed path from one vertex to another, represented by the method's two arguments. 
    inline bool path_connected(int u,int v) const {return (distance(u,v) > -1);};
    /// This method computes the set of vertices of the directed graph which have no outgoing edges and at least one incoming edge, stored in the method's argument, and returns the number of sinks found.
    unsigned int compute_sinks(std::set<int>&) const;
    /// This method computes the set of vertices of the directed graph which have no incoming edges and at least one outgoing edge, stored in the method's argument, and returns the number of sources found.
    unsigned int compute_sources(std::set<int>&) const;
    /// This method writes the directed graph out to a disk file in the DOT format for visualization by GraphViz; the method's argument is the filename.
    void write2disk(const std::string&) const override;
  };

  void Directed_Graph::compute_directedness() 
  {
    // A method to calculate how many of the edges are directed...
    unsigned int null = 0;
    std::vector<Edge>::const_iterator it;

    for(it=edges.begin(); it!=edges.end(); ++it) {
      if (it->get_direction() == Relation::disparate) null++;
    }
    number_directed = size() - null;
  }
}
#endif
