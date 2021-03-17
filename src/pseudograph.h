#include "edge.h"

#ifndef _pgraphh
#define _pgraphh

namespace SYNARMOSMA {
  class Pseudograph {
   private:
    int nvertex = 0;
    std::vector<std::pair<int,std::set<int> > > entourage;
    std::vector<Edge> edges;

    int DFS_bridge(int,int,int,int*,int*,hash_map&) const;
    int compute_bridges(hash_map&) const;
    void clear();
   public:
    Pseudograph();
    Pseudograph(int);
    ~Pseudograph();
    Pseudograph(const Pseudograph&);
    Pseudograph& operator =(const Pseudograph&);
    bool consistent() const;
    int serialize(std::ofstream&) const;
    int deserialize(std::ifstream&);
    int order() const;
    int get_candidates(std::vector<int>&) const;
    int self_loops(int) const;
    int self_loops() const;
    int degree(int) const;
    int in_degree(int) const;
    int out_degree(int) const;
    int multi_degree(int,int) const;
    int multi_degree(int,int,Relation) const;
    void add_edge(int,int,Relation);
    bool contract(int,int,Pseudograph*) const;
    bool remove(int,int,Pseudograph*) const;
    void write2disk(const std::string&,const std::vector<std::string>& = std::vector<std::string>()) const;
  };

  inline int Pseudograph::order() const
  {
    return nvertex;
  }

  inline int Pseudograph::degree(int v) const
  {
    int output = entourage[v].first + (signed) entourage[v].second.size();
    return output;
  }

  inline int Pseudograph::self_loops(int v) const
  {
    return entourage[v].first;
  }

  inline int Pseudograph::self_loops() const
  {
    int i,sum = 0;
    for(i=0; i<nvertex; ++i) {
     sum += entourage[i].first;
    }
    return sum;
  }
}
#endif

