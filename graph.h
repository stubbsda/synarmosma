#include "schema.h"
#include "edge.h"
#include "matrix.h"
#include "binary_matrix.h"
#include "pseudograph.h"

#ifndef _graphh
#define _graphh

namespace SYNARMOSMA {
  /// A class representing a graph, i.e. a one-dimensional simplicial complex, derived from the Schema class and using the Edge class for the 1-simplices.
  class Graph: public Schema {
   protected:
    /// This vector contains the list of edges of the graph, stored using the 
    /// Edge class.
    std::vector<Edge> edges;
    /// This property is a hash table that permits quick look-ups of the index 
    /// of an edge in the vector Graph::edges on the basis of the set of its two 
    /// vertices.
    hash_map index_table;

    int DFS_bridge(int,int,int,int*,int*) const;
    int DFS_cycle(int,int,std::vector<int>&,bool*) const;
    /// This method computes the eigenvalues of the graph's adjacency matrix and writes them to the argument, in ascending order.
    void adjacency_eigenvalues(std::vector<double>&) const;
    void defoliate(const Pseudograph*,std::vector<Monomial<int> >&) const;
  public:
    /// The default constructor which does nothing. 
    Graph();
    /// The standard constructor which sets the number of vertices to the argument while the number of edges is left at zero.
    Graph(int);
    Graph(std::string&);
    Graph(int,std::string&);
    Graph(int,int);
    Graph(int,double);
    Graph(const Graph&);
    Graph& operator =(const Graph&);
    /// The destructor which does nothing.
    virtual ~Graph() override;
    /// This method carries out all the operations of the Schema::clear() method and then clears the Graph::edges and Graph::index_table properties.
    virtual void clear() override;
    /// This method returns the topological energy of this graph, based on its genus.
    virtual double compute_energy() const;
    /// This method adds an edge to the graph between the two vertices specified in its two arguments, assuming an edge length of zero; the method returns true if this a new edge, false otherwise.
    virtual inline bool add_edge(int u,int v) override {return add_edge(u,v,0.0);};
    /// This method adds an edge to the graph between the two vertices specified in its first two arguments, while the third argument is the length of this edge; the method returns true if this is a new edge, false otherwise.
    virtual bool add_edge(int,int,double);
    /// This method deletes the edge connecting the two vertices specified in its arguments, returning true if there was an edge to delete and false otherwise.
    virtual bool drop_edge(int,int) override;
    /// This method deletes a vertex specified by the argument, which also requires deleting all of this vertex's edges as well; it returns true if this is successful. 
    virtual bool drop_vertex(int);
    /// This method minimizes the graph topology according to a fitness function using simulated annealing. 
    int minimize_topology(int,double,std::vector<double>&);
    /// This method renders the graph topology complete, i.e. every vertex is directly connected to every other vertex; it returns the number of edges that had to be added to the graph to achieve this result.
    int make_complete();
    /// This method first checks if the graph is consistent as an instance of the Schema class, then performs a series of sanity checks on the elements of the Graph::edges vector and returns true if everything is consistent.
    virtual bool consistent() const override;
    /// Given a set of vertices (the first argument), this method calculates the set of edges (the second argument) connecting them to the rest of the graph, i.e. the "surface" which encloses the volume (the set of vertices).
    void compute_surface(const std::set<int>&,std::set<int>&) const;
    // Hyphantic operators
    bool stellar_addition(int);
    bool stellar_deletion(int);
    virtual bool fusion(int,int);
    virtual bool foliation_x(int,int);
    virtual bool foliation_m(int,int);
    virtual int fission_x(int);
    virtual int fission_m(int);
    /// This method computes the maximum network flow from a source vertex (first argument) to a sink vertex (second argument), returning the value of the flow.
    virtual double compute_flow(int,int);
    void core(Graph*,int) const;
    int eccentricity(int) const;
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
    void vertex_centrality(std::vector<double>&,double,bool = false) const;
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
    /// This method returns the size of the graph, i.e. the number of edges.
    inline int size() const {return (signed) edges.size();};
    /// This method returns the order of the graph, i.e. the number of vertices.
    inline int order() const {return nvertex;};
    /// This method writes the graph's properties to a binary disk file and returns the number of bytes written to the file.
    virtual int serialize(std::ofstream&) const override;
    /// This method calls the clear() method on the instance and then reads the properties from a binary disk file and returns the number of bytes read.
    virtual int deserialize(std::ifstream&) override;
    /// This method writes the graph out to a disk file in the DOT format for visualization by GraphViz; the method's argument is the filename.
    virtual void write2disk(const std::string&) const;
  };
}
#endif
