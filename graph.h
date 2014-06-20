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
  double cosine_similarity(int,int) const;
  double inverse_girth() const;
  int bridge_count() const;
  int depth_first_search(int,int,int,int*,int*) const;
  double cyclicity() const;
  double connectivity() const;
  int omega() const;
  void degree_distribution(bool,std::vector<double>&) const;
  double percolation(bool) const;
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
  double return_probability(int,int) const;
  void random_walk(double*,double*,int) const;
  void compute_adjacency_matrix(Binary_Matrix*) const;
  void compute_laplacian(Matrix<double>*) const;
  void genus(int*) const;
  int minimize_topology(int,double,std::vector<double>&);
  inline int size() const {return nedge;};
};
#endif
