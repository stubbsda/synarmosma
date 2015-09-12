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

    // Returns the topological energy of this graph
    virtual double compute_energy() const;
    virtual unsigned int add_vertex();
    virtual bool drop_vertex(unsigned int);
    // A basic operator for adding an edge
    virtual bool add_edge(unsigned int,unsigned int);
    // A basic operator for undoing the above edge addition
    virtual bool drop_edge(unsigned int,unsigned int);
    // A method to minimize the graph topology according 
    // to a fitness function using simulated annealing 
    int minimize_topology(unsigned int,double,std::vector<double>&);
    // A method to render the graph topology complete
    int make_complete();
   public:
    // The usual public methods for a class
    Graph();
    Graph(const Graph&);
    Graph(const char*);
    Graph(unsigned int);
    Graph(unsigned int,unsigned int);
    Graph(unsigned int,double);
    Graph& operator =(const Graph&);
    virtual ~Graph();
    virtual void clear();
    // Hyphantic operators
    virtual bool amputation(unsigned int);
    virtual bool fusion(unsigned int,unsigned int);
    virtual bool foliation_x(unsigned int,unsigned int);
    virtual bool foliation_m(unsigned int,unsigned int);
    virtual unsigned int fission_x(unsigned int);
    virtual unsigned int fission_m(unsigned int);
    // A series of const methods to calculate various graph properties
    void core(Graph*,unsigned int) const;
    bool planar() const;
    bool biconnected() const;
    double cosine_similarity(unsigned int,unsigned int) const;
    double inverse_girth() const;
    double clustering_coefficient(unsigned int) const;
    unsigned int bridge_count() const;
    unsigned int depth_first_search(int,int,int,int*,int*) const;
    double cyclicity() const;
    double connectivity() const;
    unsigned int omega() const;
    void katz_centrality(std::vector<double>&) const;
    void degree_distribution(bool,std::vector<double>&) const;
    double percolation(bool) const;
    double cyclic_resistance() const;
    double entwinement() const;
    double completeness() const;
    unsigned int cyclomatic_number() const;
    unsigned int chromatic_number() const;
    unsigned int boundary_nodes() const;
    unsigned int max_degree() const;
    unsigned int min_degree() const;
    double average_degree() const;
    double return_probability(int,int) const;
    void random_walk(double*,double*,int) const;
    void compute_adjacency_matrix(Binary_Matrix*) const;
    void compute_laplacian(Matrix<double>*) const;
    unsigned int genus(std::vector<unsigned int>&) const;
    inline unsigned int size() const;
    inline unsigned int order() const;
    virtual void serialize(std::ofstream&) const;
    virtual void deserialize(std::ifstream&);
    friend class Nexus;
  };
}
#endif
