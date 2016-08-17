#include "complex_graph.h"

using namespace SYNARMOSMA;

Complex_Graph::Complex_Graph()
{
  nvertex = 0;
}

Complex_Graph::Complex_Graph(int order)
{
  assert(order > 0);
  nvertex = order;
  neighbours = new std::vector<int>[nvertex];
}

Complex_Graph::Complex_Graph(const Complex_Graph& source)
{
  if (nvertex > 0) delete[] neighbours;
  nvertex = source.nvertex;
  neighbours = new std::vector<int>[nvertex];
  for(int i=0; i<nvertex; ++i) {
    neighbours[i] = source.neighbours[i];
  }
}
 
Complex_Graph::~Complex_Graph()
{
  if (nvertex > 0) delete[] neighbours;
}

void Complex_Graph::add_edge(int u,int v)
{
  assert(u >= 0 && u < nvertex);
  assert(v >= 0 && v < nvertex);
  neighbours[u].push_back(v);
  neighbours[v].push_back(u);
}

int Complex_Graph::compute_loops() const
{
  int i,sum = 0;
  for(i=0; i<nvertex; ++i) {
   sum += std::count(neighbours[i].begin(),neighbours[i].end(),i); 
  }
  return sum;
}

