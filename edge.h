#include "global.h"

#ifndef _edgeh
#define _edgeh

namespace SYNARMOSMA {
  class Edge {
   protected:
    int low = -1;
    int high = -1;
    bool cyclic = false;
    double length = 0.0;
    double flow = 0.0;
    double capacity = 0.0;
    Relation direction = Relation::disparate;

    void clear();
    inline Relation get_direction(int,int) const;
   public:
    Edge();
    Edge(int,int,double = 0.0,Relation = Relation::disparate);
    Edge(const Edge&);
    virtual ~Edge();
    Edge& operator =(const Edge&);
    virtual int serialize(std::ofstream&) const;
    virtual int deserialize(std::ifstream&);
    inline void get_vertices(int*) const;
    inline void set_vertices(int,int);
    inline bool invert();
    inline int get_parity() const;
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

  int Edge::get_parity() const
  {
    if (direction == Relation::disparate) return 0;
    int output = (direction == Relation::before) ? 1 : -1;
    return output;
  }

  bool Edge::invert()
  {
    if (direction == Relation::disparate) return false;
    if (direction == Relation::before) {
      direction = Relation::after;
    }
    else {
      direction = Relation::before;
    }
    return true;
  }

  Relation Edge::get_direction(int u,int v) const
  {
#ifdef DEBUG
    assert(u == low || u == high);
    assert(v == low || v == high);
#endif
    if (direction == Relation::disparate) return direction;
    if (u < v) {
      return direction;
    }
    else {
      if (direction == Relation::before) {
        return Relation::after;
      }
      else {
        return Relation::before;
      }
    }
  }
}
#endif

