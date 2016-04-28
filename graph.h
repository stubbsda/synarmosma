/*
  Copyright 2014 Daniel Stubbs

  This file is part of Synarmosma.

  Synarmosma is free software: you can redistribute it and/or modify 
  it under the terms of the GNU General Public License as published by 
  the Free Software Foundation, either version 3 of the License, or 
  (at your option) any later version.

  Synarmosma is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with Synarmosma.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "schema.h"
#include "edge.h"
#include "matrix.h"
#include "binary_matrix.h"

#ifndef _graphh
#define _graphh

namespace SYNARMOSMA {
  // The class Graph itself
  class Graph : public Schema {
   protected:
    // The edges
    std::vector<Edge> edges;
    hash_map index_table;

    int DFS_bridge(int,int,int,int*,int*) const;
    int DFS_cycle(int,int,std::vector<int>&,bool*) const;
    void adjacency_eigenvalues(std::vector<double>&) const;
  public:
    // The usual public methods for a class
    Graph();
    Graph(int);
    Graph(int,std::string&);
    Graph(int,int);
    Graph(int,double);
    Graph(const Graph&);
    Graph& operator =(const Graph&);
    virtual ~Graph();
    virtual void clear();
    // Returns the topological energy of this graph
    virtual double compute_energy() const;
    // A basic operator for adding an edge
    bool add_edge(int,int,double = 0.0);
    // A basic operator for undoing the above edge addition
    virtual bool drop_edge(int,int);
    // A method to handle dropping a vertex, a rather complicated
    // operation for this class
    virtual bool drop_vertex(int);
    // A method to minimize the graph topology according 
    // to a fitness function using simulated annealing 
    int minimize_topology(int,double,std::vector<double>&);
    // A method to render the graph topology complete
    int make_complete();
    virtual bool consistent() const;
    // Hyphantic operators
    bool stellar_addition(int);
    bool stellar_deletion(int);
    virtual bool fusion(int,int);
    virtual bool foliation_x(int,int);
    virtual bool foliation_m(int,int);
    virtual int fission_x(int);
    virtual int fission_m(int);
    // A method to compute the maximum network flow from a source vertex to a sink vertex
    virtual double compute_flow(int,int);
    // A series of const methods to calculate various graph properties
    void core(Graph*,int) const;
    int eccentricity(int) const;
    void vertex_centrality(std::vector<double>&,double) const;
    bool planar() const;
    bool biconnected() const;
    double cosine_similarity(int,int) const;
    int girth() const;
    double inverse_girth() const;
    double clustering_coefficient(int) const;
    int bridge_count() const;
    bool bipartite() const;
    double cyclicity() const;
    double connectivity() const;
    int omega() const;
    void katz_centrality(std::vector<double>&) const;
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
    double return_probability(int,int) const;
    void random_walk(double*,double*,int) const;
    void compute_adjacency_matrix(Binary_Matrix*) const;
    void compute_laplacian(Matrix<double>*) const;
    int genus(std::vector<int>&) const;
    inline int size() const {return (signed) edges.size();};
    inline int order() const {return nvertex;};
    virtual void serialize(std::ofstream&) const;
    virtual void deserialize(std::ifstream&);
    friend class Nexus;
  };
}
#endif
