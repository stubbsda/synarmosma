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

    /// This method uses a recursive depth-first search to enumerate the number of bridges in the graph, which is the return value.
    int DFS_bridge(int,int,int,int*,int*) const;
    /// This method uses a recursive depth-first search to enumerate the number of cycles in the graph, which is the return value.
    int DFS_cycle(int,int,std::vector<int>&,bool*) const;
    /// This method computes the eigenvalues of the graph's adjacency matrix and writes them to the argument, in ascending order.
    void adjacency_eigenvalues(std::vector<double>&) const;
    /// This method is used to calculate the Tutte polynomial (the second argument) of the graph, based on the pseudograph (the first argument) which is constructed from the graph.
    void defoliate(const Pseudograph*,std::vector<Monomial<int> >&) const;
    /// This method initializes the graph to the properties of its argument.
    void initialize(const Graph&);
    /// This method accepts the number of vertices per dimension \f$n > 1\f$ (the first argument) and a vector of integers whose length is the number of dimensions \f$d > 1\f$; it returns the number of edges in the resulting lattice. Each element of the vector must be either zero (indicating a linear boundary topology) or one (toroidal boundary topology). The method builds a graph of order \f$n^d\f$ with a \f$d\f$-dimensional rectangular lattice topology, possessing a linear or toroidal boundary at each dimension according to the second argument. 
    int build_lattice(int,const std::vector<int>&);
  public:
    /// The default constructor which does nothing. 
    Graph();
    /// The standard constructor which sets the number of vertices \f$n > 0\f$ to the argument while the number of edges is left at zero.
    Graph(int);
    /// This is the constructor for "named" graphs - DÜRER, GOLOMB, HERSCHEL, PETERSEN and WAGNER - which have a fixed number of vertices and edges as well as topology.
    Graph(const std::string&);
    /// This is the contructor for a category of named graphs that also require the number of vertices \f$n > 0\f$ to be specified: COMPLETE, CHAIN, CYCLIC and CONNECTED. The latter constructs a graph on \f$n\f$ vertices by adding edges randomly until the graph is connected. 
    Graph(int,const std::string&);
    /// This constructor accepts the number of vertices \f$n > 0\f$ (the first argument) and a minimum degree \f$d > 0\f$ (the second argument) to construct a scale-free graph on \f$n\f$ vertices with minimum degree \f$d\f$.
    Graph(int,int);
    /// This constructor builds a graph with a regular lattice topology by calling the build_lattice() method using its two arguments.
    Graph(int,const std::vector<int>&);
    /// This constructor accepts the number of vertices \f$n > 0\f$ (the first argument) and a percentage \f$ 0 < \rho < 1\f$ (the second argument); this percentage is relative to the number of edges for a complete graph on \f$n\f$ vertices. 
    Graph(int,double);
    /// The standard copy constructor - it copies over all properties from the source instance to this one.
    Graph(const Graph&);
    /// The overloaded assignment operator for this class, which behaves exactly like the copy constructor for this class.
    Graph& operator =(const Graph&);
    /// The destructor which does nothing.
    virtual ~Graph() override;
    /// This method carries out all the operations of the Schema::clear() method and then clears the Graph::edges and Graph::index_table properties.
    virtual void clear() override;
    /// This method returns the topological energy of this graph, based on its genus.
    virtual double compute_energy() const;
    /// This method adds an edge to the graph between the two vertices specified in its two arguments, assuming an edge length of zero; the method returns true if this a new edge, false otherwise.
    virtual bool add_edge(int,int) override;
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
    /// Given a set of vertices (the first argument), this method calculates the set of edges (the second argument) connecting them to the rest of the graph, i.e. the "surface" which encloses the volume (the set of vertices). The method returns the number of nodes in the first argument which lie on this surface, i.e. have at least one neighbour which lies outside the set of vertices that is the first argument.
    int compute_surface(const std::set<int>&,std::set<int>&) const;
    /// This method determines if the vertex represented by its argument is part of a 3-cycle; if so, this 3-cycle is replaced by a Y topology among the three vertices by the addition of a new vertex. The method returns true if the transformation is successful, false otherwise.
    bool stellar_addition(int);
    /// This method carries out the inverse transformation of stellar_addition() - if the vertex given by the method's argument has a degree equal to three, then it is deleted and its three neighbours are placed in a 3-cycle among themselves. The method returns true if this is successful, false otherwise.
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
    /// This method computes the \f$k\f$-core of the graph and writes this new graph to its first argument, where \f$k\f$ is the second argument. The \f$k\f$-core of a graph is the set of vertices (and their accompanying edges) whose degree is greater than or equal to \f$k\f$.
    void core(Graph*,int) const;
    /// This method computes the eccentricity of a vertex, its unique argument, and returns this value; the eccentricity of a vertex is the maximum combinatorial distance between it and every other vertex in the graph.
    int eccentricity(int) const;
    /// This method returns the radius of the graph, i.e. the minimum of the vertex eccentricity over the entire set of vertices; if the graph is disconnected the method returns -1. 
    int radius() const;
    /// This method returns the diameter of the graph, i.e. the maximum of the vertex eccentricity over the entire set of vertices; if the graph is disconnected the method returns -1.
    int diameter() const;
    /// This method returns true if the graph is planar, false otherwise, using the Boost library (the Boyer-Myrvold test) to do the planarity testing.
    bool planar() const;
    /// This method returns false if there is at least one vertex whose removal disconnects the graph, otherwise it returns true.
    bool biconnected() const;
    /// This method returns the entropy of the graph, defined as the logarithm of the total number of paths of length \f$N > 0\f$, the method's argument, in the graph. This value is calculated by computing the \f$N\f$-th power of the adjacency matrix and then obtaining the sum of all its upper diagonal elements. When the argument is zero, the method uses as the length the value of Graph::nvertex less one.
    double entropy(int = 0) const;
    /// This method returns the cosine similarity between two vertices \f$u\f$ and \f$v\f$. This is defined to be the \f$(u,v)\f$ element of the square of the adjacency matrix, divided by \f$\sqrt{\deg(u)\deg(v)}\f$.
    double cosine_similarity(int,int) const;
    /// This method returns the length of the graph's largest cycle; it returns -1 if the graph is acyclic. 
    int girth() const;
    /// This method returns the inverse girth; if the graph is acyclic it returns zero.
    double inverse_girth() const;
    /// This method computes the complement of this graph, i.e. the graph with the same vertices but where two vertices are connected if and only if they are not connected in the original graph.
    void complement(Graph*) const;
    /// This method builds a pseudograph representation of the graph and then uses the defoliate() method to encode this pseudograph into a multivariate polynomial, the Tutte polynomial of the graph.
    void tutte_polynomial(std::vector<Monomial<int> >&) const;
    /// This method accepts a vertex as its argument v and computes the percentage of the distinct pairs of neighbours of v whose members are also directly connected together.
    double clustering_coefficient(int) const;
    /// This method computes the sum of the clustering coefficient over all vertices in the graph and divides this by the number of vertices, which is then returned.
    double clustering_coefficient() const;
    /// This method sums the combinatorial distance between every distinct pair of vertices in the graph and divides this by the number of such pairs, returning the result.
    double mean_path_length() const;
    /// This method returns the number of edges which are bridges, i.e. edges which if deleted would disconnect the graph.
    int bridge_count() const;
    /// This method returns false if the graph contains a cycle of odd length and true otherwise.
    bool bipartite() const;
    /// This method returns the number of graph vertices which are \f$N\ge 0\f$ steps or less from the vertex which is the method's first argument, where \f$N\f$ is the second argument. 
    int compactness(int,int) const;
    /// This method returns the graph's cyclicity, i.e. the percentage of edges whose deletion would not disconnect the graph.
    double cyclicity() const;
    /// This method returns the value of the graph connectivity, defined as the sum over all edges of \f$1/\sqrt{\deg(u)\deg(v)}\f$, where \f$u\f$ and \f$v\f$ are the two vertices joined by the edge. 
    double connectivity() const;
    /// This method returns the value of the algebraic connectivity of the graph, i.e. the smallest non-zero eigenvalue of graph Laplacian.
    double algebraic_connectivity() const;
    /// This method computes the eigenvector centrality, a vector whose length is equal to the graph's order and which measures each vertex's centrality or importance in the graph. The output is the method's first argument, while the second argument sets the tolerance for the recursion process that generates the centrality vector and the final argument determines whether or not Katz centrality is calculated.
    void vertex_centrality(std::vector<double>&,double,bool = false) const;
    /// This method calculates a histogram of vertex degrees which is stored in the method's second argument. The first argument controls whether the binning is logarithmic or not; if true the bin width doubles, 1, 2, 4, 8, 16 and so forth.
    void degree_distribution(bool,std::vector<double>&) const;
    /// This method computes the percentage of randomly selected vertices (site percolation) or edges (bond percolation) that need to be removed from the graph to eliminate its "giant component" and returns this percentage. If the argument is true the method uses site percolation, otherwise bond percolation. 
    double percolation(bool) const;
    /// This method returns the index of graph complexity discussed in the paper Klein et al., "Graph Cyclicity, Excess Conductance and Resistance Deficit", J. Math. Chem., 30:271-287 (2001).
    double cyclic_resistance() const;
    /// This method returns the median degree of the graph, i.e. the degree \f$d\f$ for which half the vertices in the graph have a degree less than or equal to \f$d\f$. This is in general a floating point number to account for the fact that it is unlikely that the sum of degree histogram bins is exactly half the vertices in the graph so the median has to be corrected. 
    double median_degree() const;
    /// This method returns the graph's entwinement. This is defined to be \f$(\lambda_\textrm{max} - \kappa)/(d_\textrm{max} - \kappa)\f$, where \f$\lambda_\textrm{max}\f$ is the largest eigenvalue of its adjacency matrix, \f$d_\textrm{max}\f$ the maximum degree and \f$\kappa = \max(d_\textrm{avg},\sqrt{d_\textrm{max}})\f$. 
    double entwinement() const;
    /// This method returns the graph's completeness, i.e. the graph's size divided by \f$N(N-1)/2\f$, where \f$N\f$ is the graph's order; this measures how closely it approximates the complete graph on \f$N\f$ vertices.
    double completeness() const;
    /// This method returns the graph's circuit rank, defined to be the number of graph components minus the number of vertices plus the number of edges.
    int circuit_rank() const;
    /// This method returns the graph's chromatic number, i.e. the smallest number of colours which can be used to colour the vertices such that no two neighbouring vertices share the same colour.
    int chromatic_number() const;
    /// This method returns the minimum degree of the graph.
    int max_degree() const;
    /// This method returns the maximum degree of the graph.
    int min_degree() const;
    /// This method returns the arithmetic mean of the vertex degrees of the graph.
    double average_degree() const;
    /// This method returns true if the graph is Eulerian, i.e. it is connected and all of its vertices have even degree, and false otherwise.
    virtual bool eulerian() const;
    /// This method attempts to construct a Hamiltonian path in the graph (which must be connected) and which, if successful, is stored as a vector of vertices in the third argument. If the first argument is true the method attempts to construct a Hamiltonian cycle, while the second argument contains the number of attempts to be made. The optional fourth argument is the initial vertex; if it isn't provided a vertex is chosen at random. The method returns the number of attempts which were made to construct the path. 
    virtual int compute_hamiltonian_path(bool,int,std::vector<int>&,int = -1) const;
    /// This method's argument represent a base vertex v, a length l and the number of trials n - the method carries out n random walks from v of length l and returns the percentage of these walks which return at least once to the starting vertex v.
    double return_probability(int,int,int = 100) const;
    /// This method uses the return_probability() method to do a general analysis of random walks on the graph topology. Its arguments are the length of the random walk, the percentage of vertices randomly selected that are used as a starting point and the number of walks per starting point. The method returns a pair of doubles, the average return probability and its standard deviation over the number of starting points.
    std::pair<double,double> random_walk(int,double,int = 15) const;
    /// This method computes the adjacency matrix of the graph and stores it in an instance of the Binary_Matrix class, which is the method's unique argument.
    void compute_adjacency_matrix(Binary_Matrix*) const;
    /// This method computes the deformed graph Laplacian, defined as \f$\Delta_G(s) = I + s^2 (D - I) - s A\f$, where \f$I\f$ is the identity matrix, \f$D\f$ is the diagonal matrix of vertex degrees, \f$A\f$ is the adjacency matrix and \f$s\f$ is a complex number. It is easy to see that \f$\Delta_G(1) = L_G\f$, the standard graph Laplacian, while \f$\Delta_G(0) = I\f$. The method returns the number of non-zero entries in this deformed Laplacian.
    int compute_deformed_laplacian(std::complex<double>,Matrix<std::complex<double> >*) const;
    /// This method computes the graph Laplacian, defined as \f$ L_G = D - A\f$ where \f$D\f$ is the diagonal matrix of vertex degrees and \f$A\f$ is the adjacency matrix of the graph; the method returns the number of non-zero entries in the graph Laplacian.
    int compute_laplacian(Matrix<double>*) const;
    /// This method returns the value of the graph's genus; if the graph is planar this is just zero. For a non-planar graph, the method computes an upper and lower bound for the genus; if these coincide the method returns this value, otherwise it puts the two values in the method's argument and returns -1.
    int genus(std::vector<int>&) const;
    /// This method returns the size of the graph, i.e. the number of edges.
    int size() const;
    /// This method writes the graph's properties to a binary disk file and returns the number of bytes written to the file.
    virtual int serialize(std::ofstream&) const override;
    /// This method calls the clear() method on the instance and then reads the properties from a binary disk file and returns the number of bytes read.
    virtual int deserialize(std::ifstream&) override;
    /// This method writes the graph out to a disk file in the DOT format for visualization by GraphViz; the method's first argument is the filename, the optional second argument is a vector of strings which are labels for the graph's vertices.
    virtual void write2disk(const std::string&,const std::vector<std::string>& = std::vector<std::string>()) const;
    /// This operator returns the tensor product of its two arguments; this product is a graph whose order is the square of the operands' order (which must be the same) and which has an edge between two vertices when the two pairs of vertices \f$(g_1,h_1)\f$ and \f$(g_2,h_2)\f$ are both neighbours in the respective operand graphs \f$G\f$ and \f$H\f$. Hence there is an edge in \f$G\otimes H\f$ when \f$g_1\f$ is a neighbour of \f$g_2\f$ in \f$G\f$ and \f$h_1\f$ is a neighbour of \f$h_2\f$ in \f$H\f$. 
    friend Graph operator *(const Graph&,const Graph&);
  };

  inline int Graph::size() const {
    return (signed) edges.size(); 
  }

  inline double Graph::inverse_girth() const
  {
    int output = girth();
    double g = 0.0;
    if (output > 0) g = 1.0/double(output);
  
    return g;
  }

  inline int Graph::circuit_rank() const
  {
    int chi = size() - nvertex;
    if (connected()) return 1 + chi;
    std::vector<int> components;
    int n = component_analysis(components);
    return n + chi;
  }

  inline int Graph::max_degree() const
  {
    int i,n,output = 0;
    for(i=0; i<nvertex; ++i) {
      n = (signed) neighbours[i].size();
      if (n > output) output = n;
    }
    return output;
  }

  inline int Graph::min_degree() const
  {
    int i,n,output = 1 + size();
    for(i=0; i<nvertex; ++i) {
      n = (signed) neighbours[i].size();
      if (n < output) output = n;
    }
    return output;
  }

  inline double Graph::average_degree() const
  {
    int i;
    double sum = 0.0;
    for(i=0; i<nvertex; ++i) {
      sum += double(neighbours[i].size());
    }
    return sum/double(nvertex);
  }

  inline bool Graph::eulerian() const
  {
    if (!connected()) return false;

    for(int i=0; i<nvertex; ++i) {
      if (neighbours[i].size() % 2 != 0) return false;
    }
    return true;
  }

  inline double Graph::completeness() const
  {
    double output = 0.0;
    if (nvertex == 1) return output;

    output = 2.0*double(size())/double(nvertex*(nvertex-1));
    return output;
  }
}
#endif
