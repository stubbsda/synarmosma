#include "edge.h"

#ifndef _pgraphh
#define _pgraphh

namespace SYNARMOSMA {
  /// A class representing a pseudograph, i.e. a graph which permits multi-edges and self-loops.
  class Pseudograph {
   private:
    /// The number of vertices in the pseudograph, which is also the 
    /// size of the array Pseudograph::entourage that stores the topology  
    /// of the pseudograph.
    int nvertex = 0;
    /// A vector of length equal to Pseudograph::nvertex, with 
    /// each element a pair consisting of an integer, the number of 
    /// self-loops for this vertex, and a set of integers, storing 
    /// the index of elements of Pseudograph::edges which join this 
    /// vertex to another.  
    std::vector<std::pair<int,std::set<int> > > entourage;
    /// This vector contains the list of edges of the pseudograph, stored using the 
    /// Edge class; there may be multiple edges with the same direction joining together 
    /// two vertices since this is a pseudograph. 
    std::vector<Edge> edges;

    /// This method performs a depth-first search for bridges in the pseudograph - the first two arguments are the vertices, the third is the depth, the next two integer arrays store which vertices have already been seen and finally a hash map to index the bridges found; the return value is the number of bridges found.
    int DFS_bridge(int,int,int,int*,int*,hash_map&) const;
    /// This method is the one normally called to carry out a complete analysis of bridges in the pseudograph and uses the DFS_bridge() method after initializing the arrays and the hash map; the return value is the total number of bridges in the pseudograph.
    int compute_bridges(hash_map&) const;
    /// This method sets Pseudograph::nvertex to zero and empties the Pseudograph::entourage and Pseudograph::edges vectors. 
    void clear();
   public:
    /// The default constructor which does nothing.
    Pseudograph();
    /// This constructor takes the argument as the pseudograph's order, which no edges or self-loops at all. 
    Pseudograph(int);
    /// The destructor for this class, which does nothing.
    ~Pseudograph();
    /// The standard copy constructor that copies over the value of the Pseudograph::nvertex, Pseudograph::entourage and Pseudograph::edges properties.
    Pseudograph(const Pseudograph&);
    /// The overloaded assignment operator that behaves exactly like the copy constructor for this class.
    Pseudograph& operator =(const Pseudograph&);
    /// This method verifies that all of the values of the Pseudograph::entourage and Pseudograph::edges properties make sense, returning true if this is so and false otherwise. 
    bool consistent() const;
    /// This method writes the instance properties to a binary disk file and returns the number of bytes written to the file.
    int serialize(std::ofstream&) const;
    /// This method calls the clear() method on the instance and then reads the properties from a binary disk file and returns the number of bytes read.
    int deserialize(std::ifstream&);
    /// This method simply returns the value of Pseudograph::nvertex, the order of the pseudograph.
    int order() const;
    /// This method uses the compute_bridges() method to compute all the bridges in the pseudograph - the output value of the method - and so assemble a list of candidate edges for contraction, i.e. those which are not bridges.
    int get_candidates(std::vector<int>&) const;
    /// This method returns the number of self-loops for a given vertex, the method's argument n, so the value of Pseudograph::entourages[n].first.
    int self_loops(int) const;
    /// This method calculates and returns the number of self-loops in this pseudograph.
    int self_loops() const;
    /// This method returns the total number of edges connected to the vertex which is the method's argument.
    int degree(int) const;
    /// This method returns the in degree of the vertex which is its argument, i.e. the number of directed edges which terminate at this vertex.
    int in_degree(int) const;
    /// This method returns the out degree of the vertex which is its argument, i.e. the number of directed edges which originate at this vertex.
    int out_degree(int) const;
    /// This method returns the total number of edges joining together the two vertices given by the method's arguments. 
    int multi_degree(int,int) const;
    /// This method returns the total number of edges of a given direction joining together the two vertices given by the method's arguments. 
    int multi_degree(int,int,Relation) const;
    /// This method adds an edge joining the two vertices which are its first arguments, with diection given by the final argument, to the pseudograph.
    void add_edge(int,int,Relation);
    /// This method fuses together two vertices (the first two arguments), removing any edges between them; it functions by replacing the greater of the two vertices by the lesser. The final argument is the output pseudograph in which the vertices have been fused.
    bool contract(int,int,Pseudograph*) const;
    /// This method removes any edges between the two vertices that are the first arguments of the method, writing the new output pseudograph as the final argument.
    bool remove(int,int,Pseudograph*) const;
    /// This method writes the pseudograph to a disk file in the DOT format for visualization by GraphViz; the method's first argument is the filename, the optional second argument is a vector of strings which are labels for the pseudograph's vertices.
    void write2disk(const std::string&,const std::vector<std::string>& = std::vector<std::string>()) const;
  };

  inline int Pseudograph::order() const
  {
    return nvertex;
  }

  inline int Pseudograph::degree(int v) const
  {
    if (v < 0 || v >= nvertex) throw std::invalid_argument("The Pseudograph::degree argument must be a vertex of the pseudograph!");
    int output = entourage[v].first + (signed) entourage[v].second.size();
    return output;
  }

  inline int Pseudograph::self_loops(int v) const
  {
    if (v < 0 || v >= nvertex) throw std::invalid_argument("The Pseudograph::self_loops argument must be a vertex of the pseudograph!");
    return entourage[v].first;
  }

  inline int Pseudograph::self_loops() const
  {
    int i,sum = 0;
    for(i=0; i<nvertex; ++i) {
     sum += entourage[i].first;
    }
    return sum;
  }
}
#endif

