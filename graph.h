#include "schema.h"
#include "matrix.h"
#include "binary_matrix.h"

#ifndef _graphh
#define _graphh

// The class Graph itself
class Graph : public Schema {
 protected:
  // The number of edges
  int nedge;
  // Returns the topological energy of this graph
  virtual double compute_energy() const;

 public:
  Graph();
  Graph(const char*);
  Graph(int);
  Graph(int,double);
  virtual ~Graph();
  virtual void clear();
  bool planar() const;
  bool biconnected() const;
  double inverse_girth() const;
  int bridge_count() const;
  int depth_first_search(int,int,int,int*,int*) const;
  double cyclicity() const;
  double connectivity() const;
  int omega() const;
  double return_probability(int,int) const;
  double cyclic_resistance() const;
  double entwinement() const;
  double completeness() const;
  int cyclomatic_number() const;
  int chromatic_number() const;
  int boundary_nodes() const;
  int max_degree() const;
  int min_degree() const;
  double average_degree() const;
  bool add_edge(int,int);
  void compute_adjacency_matrix(Binary_Matrix*) const;
  void genus(int*) const;
  int minimize_topology(int,double,std::vector<double>&);
  void build_laplacian(Matrix<double>*);
  inline int size() const {return nedge;};
  friend void build_graph(int,int,Graph*);
};
#endif
