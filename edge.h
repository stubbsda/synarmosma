#include "global.h"

#ifndef _edgeh
#define _edgeh

namespace SYNARMOSMA {
  class Edge {
   protected:
    double length;
    bool cyclic;
    double flow;
    double capacity;
    RELATION direction;
    std::set<int> vertices;

    void clear();
    void serialize(std::ofstream&) const;
    void deserialize(std::ifstream&);
   public:
    Edge();
    Edge(int,int,double = 0.0,RELATION = DISPARATE);
    Edge(const Edge&);
    ~Edge();
    Edge& operator =(const Edge&);
    inline void get_vertices(int*) const;
    inline void set_vertices(int u,int v);
    friend class Graph;
    friend class Directed_Graph;
    friend class Logic_Graph;
  };

  void Edge::set_vertices(int u,int v) 
  {
#ifdef DEBUG
    assert(u != v);
#endif
    vertices.clear();
    vertices.insert(u);
    vertices.insert(v);
  }

  void Edge::get_vertices(int* vx) const
  {
    int i = 0;
    std::set<int>::const_iterator it;
    for(it=vertices.begin(); it!=vertices.end(); ++it) {
      vx[i] = *it; ++i;
    }
  }
}
#endif

