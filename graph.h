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
    /// This method fuses together the two vertices that are its arguments; if the final argument is -1, it chooses a neighbouring vertex of the first argument at random to fuse with this first argument. In this fusion the edge sets of the two vertices are combined together and the method returns true if the fusion is successful.
    virtual bool fusion(int,int);
    /// This method eliminates the edge joining the two vertices that are its arguments; if the final argument is -1, it chooses a neighbouring vertex of the first argument at random to eliminate, returning true if successful.
    virtual bool foliation_x(int,int);
    /// This method adds an edge joining the two vertices that are its arguments; if the final argument is -1, it chooses another distinct vertex at random and adds an edge between them, returning true if this is a new edge.
    virtual bool foliation_m(int,int);
    /// This method adds a new vertex to the graph and an edge joining this new vertex to the vertex that is the argument of the method, returning the index of the newly created vertex.
    virtual int fission_x(int);
    /// This method adds a new vertex to the graph along with edges joining this new vertex to the vertex that is the argument of this method and all its neighbours, while returning the index of the newly created vertex. 
    virtual int fission_m(int);
    /// This method computes the maximum network flow from a source vertex (first argument) to a sink vertex (second argument), returning the value of the flow.
    virtual double compute_flow(int,int);
    void core(Graph*,int) const;
    int eccentricity(int) const;
    /// This method returns true if the graph is planar, false otherwise. 
    bool planar() const;
    bool biconnected() const;
    double cosine_similarity(int,int) const;
    /// This method computes the length of the graph's largest cycle; it returns -1 if the graph is acyclic. 
    int girth() const;
    /// This method computes the inverse girth; if the graph is acyclic it returns zero.
    inline double inverse_girth() const;
    /// This method computes the complement of this graph, i.e. the graph with the same vertices but where two vertices are connected if and only if they are not connected in the original graph.
    void complement(Graph*) const;
    void tutte_polynomial(std::vector<Monomial<int> >&) const;
    double clustering_coefficient(int) const;
    double clustering_coefficient() const;
    double mean_path_length() const;
    /// This method computes the number of edges which are bridges, i.e. edges which if deleted would disconnect the graph.
    int bridge_count() const;
    bool bipartite() const;
    /// This method computes the graph's cyclicity, i.e. the percentage of edges whose deletion would not disconnect the graph.
    double cyclicity() const;
    double connectivity() const;
    int omega() const;
    double algebraic_connectivity() const;
    void vertex_centrality(std::vector<double>&,double,bool = false) const;
    void degree_distribution(bool,std::vector<double>&) const;
    double percolation(bool) const;
    double cyclic_resistance() const;
    double entwinement() const;
    /// This method computes the graph's completeness, i.e. the graph's size divided by \f$N(N-1)/2\f$, where \f$N\f$ is the graph's order; this measures how closely it approximates the complete graph on \f$N\f$ vertices.
    inline double completeness() const;
    /// This method computes the graph's circuit rank, defined to be the number of graph components minus the number of vertices plus the number of edges.
    inline int circuit_rank() const;
    /// This method computes the graph's chromatic number, i.e. the smallest number of colours which can be used to colour the vertices such that no two neighbouring vertices share the same colour.
    int chromatic_number() const;
    int boundary_nodes() const;
    /// This method returns the minimum degree of the graph.
    inline int max_degree() const;
    /// This method returns the maximum degree of the graph.
    inline int min_degree() const;
    /// This method returns the arithmetic mean of the vertex degrees of the graph.
    inline double average_degree() const;
    double return_probability(int,int) const;
    void random_walk(double*,double*,int) const;
    /// This method computes the adjacency matrix of the graph and stores it in an instance of the Binary_Matrix class, which is the method's unique argument.
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

  double Graph::inverse_girth() const
  {
    int output = girth();
    double g = 0.0;
    if (output > 0) g = 1.0/double(output);
  
    return g;
  }

  int Graph::circuit_rank() const
  {
    int chi = size() - nvertex;
    if (connected()) return 1 + chi;
    std::vector<int> components;
    int n = component_analysis(components);
    return n + chi;
  }

  int Graph::max_degree() const
  {
    int i,n,output = 0;
    for(i=0; i<nvertex; ++i) {
      n = (signed) neighbours[i].size();
      if (n > output) output = n;
    }
    return output;
  }

  int Graph::min_degree() const
  {
    int i,n,output = 1 + size();
    for(i=0; i<nvertex; ++i) {
      n = (signed) neighbours[i].size();
      if (n < output) output = n;
    }
    return output;
  }

  double Graph::average_degree() const
  {
    int i;
    double sum = 0.0;
    for(i=0; i<nvertex; ++i) {
      sum += double(neighbours[i].size());
    }
    return sum/double(nvertex);
  }

  double Graph::completeness() const
  {
    double output = 0.0;
    if (nvertex == 1) return output;

    output = 2.0*double(size())/double(nvertex*(nvertex-1));
    return output;
  }
}
#endif
