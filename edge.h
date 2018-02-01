#include "global.h"

#ifndef _edgeh
#define _edgeh

namespace SYNARMOSMA {
  class Edge {
   protected:
    int low;
    int high;
    double length;
    bool cyclic;
    double flow;
    double capacity;
    int direction;

    void clear();
    inline int get_direction(int,int) const;
   public:
    Edge();
    Edge(int,int,double = 0.0,int = UNDIRECTED);
    Edge(const Edge&);
    virtual ~Edge();
    Edge& operator =(const Edge&);
    int serialize(std::ofstream&) const;
    int deserialize(std::ifstream&);
    inline void get_vertices(int*) const;
    inline void set_vertices(int,int);
    friend class Graph;
    friend class Directed_Graph;
    friend class Logic_Graph;
  };

  void Edge::set_vertices(int u,int v) 
  {
#ifdef DEBUG
    assert(u != v);
#endif
    if (u < v) {
      low = u; high = v;
    }
    else {
      low = v; high = u;
    }
  }

  void Edge::get_vertices(int* vx) const
  {
    vx[0] = low; vx[1] = high;
  }

  int Edge::get_direction(int u,int v) const
  {
#ifdef DEBUG
    assert(u == low || u == high);
    assert(v == low || v == high);
#endif
    if (direction == UNDIRECTED) return direction;
    if (u < v) {
      return direction;
    }
    else {
      return -direction;
    }
  }
}
#endif

