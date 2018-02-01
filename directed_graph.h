#include "graph.h"

#ifndef __dgraph
#define __dgraph

namespace SYNARMOSMA {
  class Directed_Graph: public Graph {
   protected:
    int number_directed;
    bool directed_cycle(const std::vector<int>&,int,int) const;

   public:
    Directed_Graph();
    Directed_Graph(int);
    Directed_Graph(int,double);
    Directed_Graph(const Directed_Graph&);
    Directed_Graph& operator =(const Directed_Graph&);
    virtual ~Directed_Graph();
    virtual int serialize(std::ofstream&) const;
    virtual int deserialize(std::ifstream&);
    bool add_edge(int,int,int,double = 0.0);
    // A method to compute the maximum network flow from a source vertex 
    // to a sink vertex
    virtual double compute_flow(int,int);
    virtual int distance(int,int) const;
    virtual void compute_distances(edge_hash&) const; 
    bool mutate_edge(int,int);
    void compute_directedness();
    bool acyclic() const;
    bool path_connected(int,int) const;
    void compute_sinks(std::set<int>&) const;
    void compute_sources(std::set<int>&) const;
    friend class Propositional_System;
  };
}
#endif
