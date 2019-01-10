#include "schema.h"
#include "edge.h"
#include "matrix.h"
#include "binary_matrix.h"
#include "pseudograph.h"

#ifndef _graphh
#define _graphh

namespace SYNARMOSMA {
  // The class Graph itself
  class Graph: public Schema {
   protected:
    // The edges
    std::vector<Edge> edges;
    hash_map index_table;

    int DFS_bridge(int,int,int,int*,int*) const;
    int DFS_cycle(int,int,std::vector<int>&,bool*) const;
    void adjacency_eigenvalues(std::vector<double>&) const;
    void defoliate(const Pseudograph*,std::vector<Monomial<int> >&) const;
  public:
    // The usual public methods for a class
    Graph();
    Graph(int);
    Graph(std::string&);
    Graph(int,std::string&);
    Graph(int,int);
    Graph(int,double);
    Graph(const Graph&);
    Graph& operator =(const Graph&);
    ~Graph() override;
    void clear() override;
    // Returns the topological energy of this graph
    virtual double compute_energy() const;
    // A basic operator for adding an edge
    bool add_edge(int,int) override;
    virtual bool add_edge(int,int,double);
    // A basic operator for undoing the above edge addition
    bool drop_edge(int,int) override;
    // A method to handle dropping a vertex, a rather complicated
    // operation for this class
    virtual bool drop_vertex(int);
    // A method to minimize the graph topology according 
    // to a fitness function using simulated annealing 
    int minimize_topology(int,double,std::vector<double>&);
    // A method to render the graph topology complete
    int make_complete();
    bool consistent() const override;
    // Given a set of vertices, calculate the set of edges connecting them to the 
    // rest of the graph, i.e. the "surface" which encloses the volume (the set of vertices)
    void compute_surface(const std::set<int>&,std::set<int>&) const;
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
    void complement(Graph*) const;
    void tutte_polynomial(std::vector<Monomial<int> >&) const;
    double clustering_coefficient(int) const;
    double clustering_coefficient() const;
    double mean_path_length() const;
    int bridge_count() const;
    bool bipartite() const;
    double cyclicity() const;
    double connectivity() const;
    int omega() const;
    double algebraic_connectivity() const;
    void katz_centrality(std::vector<double>&) const;
    void degree_distribution(bool,std::vector<double>&) const;
    double percolation(bool) const;
    double cyclic_resistance() const;
    double entwinement() const;
    double completeness() const;
    int circuit_rank() const;
    int chromatic_number() const;
    int boundary_nodes() const;
    int max_degree() const;
    int min_degree() const;
    double average_degree() const;
    double return_probability(int,int) const;
    void random_walk(double*,double*,int) const;
    void compute_adjacency_matrix(Binary_Matrix*) const;
    int compute_deformed_laplacian(std::complex<double>,Matrix<std::complex<double> >*) const;
    int compute_laplacian(Matrix<double>*) const;
    int genus(std::vector<int>&) const;
    inline int size() const {return (signed) edges.size();};
    inline int order() const {return nvertex;};
    int serialize(std::ofstream&) const override;
    int deserialize(std::ifstream&) override;
    virtual void write2disk(const std::string&) const;
  };
}
#endif
