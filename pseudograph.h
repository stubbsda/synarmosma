#include "global.h"

#ifndef _cgraphh
#define _cgraphh

namespace SYNARMOSMA {
  class Pseudograph {
   private:
    int nvertex;
    hash_map bridge_index;
    std::vector<int>* neighbours;

    int DFS_bridge(int,int,int,int*,int*,hash_map&) const;
    int compute_bridges(hash_map&) const;
   public:
    Pseudograph();
    Pseudograph(int);
    ~Pseudograph();
    Pseudograph(const Pseudograph&);
    Pseudograph& operator =(const Pseudograph&);
    int get_candidates(std::vector<int>&) const;
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

