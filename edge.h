#include "global.h"

#ifndef _edgeh
#define _edgeh

class Edge {
 private:
  double length;
  std::string name;
  std::vector<std::string> entourage;
  int direction;
  bool cyclic;
  double flow;
  double capacity;
  unsigned int nodes[2];
  unsigned int colour;

  void clear();
 public:
  Edge();
  Edge(const Edge&);
  ~Edge();
  Edge& operator=(const Edge&);
};
#endif
