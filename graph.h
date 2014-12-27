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
#include "matrix.h"
#include "binary_matrix.h"

#ifndef _graphh
#define _graphh

namespace SYNARMOSMA {
  // The class Graph itself
  class Graph : public Schema {
   protected:
    // The number of edges
    int nedge;

    // Returns the topological energy of this graph
    virtual double compute_energy() const;
    // A basic operator for adding an edge
    virtual bool add_edge(int,int);
    // A method to minimize the graph topology according 
    // to a fitness function using simulated annealing 
    int minimize_topology(int,double,std::vector<double>&);
   public:
    // The usual public methods for a class
    Graph();
    Graph(const Graph&);
    Graph(const char*);
    Graph(int);
    Graph(int,int);
    Graph(int,double);
    virtual ~Graph();
    virtual void clear();
    // Hyphantic operators
    virtual bool amputation(int);
    virtual bool fusion(int,int);
    virtual bool foliation_x(int,int);
    virtual bool foliation_m(int,int);
    virtual int fission_x(int);
    virtual int fission_m(int);
    // A series of const methods to calculate various graph properties
    void core(Graph*,int) const;
    bool planar() const;
    bool biconnected() const;
    double cosine_similarity(int,int) const;
    double inverse_girth() const;
    double clustering_coefficient(int) const;
    int bridge_count() const;
    int depth_first_search(int,int,int,int*,int*) const;
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
    inline int size() const {return nedge;};
    inline int order() const {return nvertex;};
    friend class Nexus;
  };
}
#endif
