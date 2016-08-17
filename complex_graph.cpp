#include "complex_graph.h"

using namespace SYNARMOSMA;

Complex_Graph::Complex_Graph()
{
  n = 0;
}

Complex_Graph::Complex_Graph(int order)
{
  assert(order > 0);
  n = order;
  neighbours = new std::vector<int>[n];
}

Complex_Graph::~Complex_Graph()
{
  if (n > 0) delete[] neighbours;
}

void Complex_Graph::add_edge(int u,int v)
{
  assert(u >= 0 && u < n);
  assert(v >= 0 && v < n);
  neighbours[u].push_back(v);
  neighbours[v].push_back(u);
}

