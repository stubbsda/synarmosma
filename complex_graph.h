#include "global.h"

#ifndef _cgraphh
#define _cgraphh

namespace SYNARMOSMA {
  class Complex_Graph {
   private:
    int nvertex;
    hash_map bridge_index;
    std::vector<int>* neighbours;

    int DFS_bridge(int,int,int,int*,int*,hash_map&) const;
    int compute_bridges(hash_map&) const;
   public:
    Complex_Graph();
    Complex_Graph(int);
    ~Complex_Graph();
    Complex_Graph(const Complex_Graph&);
    Complex_Graph& operator =(const Complex_Graph&);
    int get_candidates(std::vector<int>&) const;
    inline int get_loops() const;
    void add_edge(int,int);
    void contract(int,int,Complex_Graph*) const;
    void remove(int,int,Complex_Graph*) const;
    void display() const;
    friend class Graph;
  };

  int Complex_Graph::get_loops() const
  {
    int i,sum = 0;
    for(i=0; i<nvertex; ++i) {
     sum += std::count(neighbours[i].begin(),neighbours[i].end(),i); 
    }
    return sum;
  }
}
#endif

