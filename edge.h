#include "global.h"

#ifndef _edgeh
#define _edgeh

enum DIRECTION 
{
    FORWARD,
    BACKWARD
};

class Edge {
 private:
  double length;
  DIRECTION arrow;
  bool cyclic;
  double flow;
  double capacity;
  int nodes[2];

  void clear();
 public:
  Edge();
  Edge(int,int);
  Edge(int,int,DIRECTION);
  Edge(const Edge&);
  ~Edge();
  Edge& operator =(const Edge&);
  inline std::string key() const {std::string k = make_key(nodes[0],nodes[1]); k += (arrow == FORWARD) ? ":1" : ":-1"; return k;};
  friend class Directed_Graph;
};
#endif
