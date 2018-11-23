#include "global.h"

#ifndef _cgraphh
#define _cgraphh

namespace SYNARMOSMA {
  /// A class representing a pseudograph, i.e. a graph which permits multi-edges and self-loops.
  class Pseudograph {
   private:
    /// The number of vertices in the pseudograph, which is 
    /// also the size of the array of vectors that stores the 
    /// neighbour lists.
    int nvertex = 0;
    std::vector<int>* neighbours;

    int DFS_bridge(int,int,int,int*,int*,hash_map&) const;
    int compute_bridges(hash_map&) const;
    void clear();
   public:
    /// The default constructor which does nothing.
    Pseudograph();
    /// A constructor which accepts the order of the pseudograph as its unique argument, setting the value of Pseudograph::nvertex and allocating the memory for Pseudograph::neighbours. 
    Pseudograph(int);
    ~Pseudograph();
    Pseudograph(const Pseudograph&);
    Pseudograph& operator =(const Pseudograph&);
    int serialize(std::ofstream&) const;
    int deserialize(std::ifstream&);
    int get_candidates(std::vector<int>&) const;
    /// This method calculates and returns the number of self-loops in this pseudograph.
    inline int get_loops() const;
    void add_edge(int,int);
    void contract(int,int,Pseudograph*) const;
    void remove(int,int,Pseudograph*) const;
    friend class Graph;
  };

  int Pseudograph::get_loops() const
  {
    int i,sum = 0;
    for(i=0; i<nvertex; ++i) {
     sum += std::count(neighbours[i].begin(),neighbours[i].end(),i); 
    }
    return sum;
  }
}
#endif

