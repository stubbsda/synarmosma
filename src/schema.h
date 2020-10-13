#include "global.h"

#ifndef _schemah
#define _schemah

namespace SYNARMOSMA {
  /// A class representing a minimal one-dimensional skeleton: a collection of vertices labelled by integers that are connected by a neighbour relation.
  class Schema {
   protected:
    /// This integer property stores the number of vertices in this 
    /// 1-skeleton and is incremented by the add_vertex() method.
    int nvertex = 0;
    /// This vector is of length equal to Schema::nvertex and contains 
    /// the set of neighbours for the vertex in question.  
    std::vector<std::set<int> > neighbours;

   public:
    /// The default constructor which does nothing.
    Schema();
    /// This constructor sets the Schema::nvertex property to its argument and adds the appropriate number of empty sets to the vector Schema::neighbours.
    Schema(int);
    /// The standard copy constructor which copies over the values of Schema::nvertex and Schema::neighbours.
    Schema(const Schema&);
    /// The standard overloaded assignment operator which copies over the values of Schema::nvertex and Schema::neighbours.
    Schema& operator =(const Schema&);
    /// The destructor which does nothing.
    virtual ~Schema();
    /// This method writes the instance properties to a binary disk file and returns the number of bytes written to the file.
    virtual int serialize(std::ofstream&) const;
    /// This method calls the clear() method on the instance and then reads the properties from a binary disk file and returns the number of bytes read.
    virtual int deserialize(std::ifstream&);
    /// This method determines if the 1-skeleton is connected (i.e. every vertex can be reached by every other vertex) and returns true if this is the case.
    virtual bool connected() const;
    /// This method checks if there is a direct connection between the two vertex arguments, i.e. if an edge connects these two vertices, returning true if so.
    virtual bool connected(int,int) const;
    /// This method restores the Schema to its default state, with Schema::nvertex equal to zero and Schema::neighbours empty.
    virtual void clear();
    /// This method adds a vertex to the Schema, incrementing Schema::nvertex and adding a new empty set to Schema::neighbours; it returns the index of the newly created vertex.
    int add_vertex();
    /// This method adds an edge between the two vertex arguments, returning false if such an edge exists already and true otherwise.
    virtual bool add_edge(int,int);
    /// This method removes the edge between the two vertex arguments, returning false if no such edge exists and true otherwise.
    virtual bool drop_edge(int,int);
    /// This method computes the combinatorial distance between the two vertex arguments, i.e. the smallest number of "hops" to get from the first argument to the second; if there is no such path (only possible if connected() is false), the method returns -1. 
    virtual int distance(int,int) const;
    /// This method computes a topological geodesic between two distinct vertices (source and target, the method's first two arguments), i.e. a path whose length is equal to the value returned by the distance() method. The path is stored in the method's final argument, as a sequence of neighbouring vertices whose final element is the target vertex (second argument); if no path exists the vector will be empty upon return. 
    virtual void compute_shortest_path(int,int,std::vector<int>&) const;
    /// This method computes the complete set of combinatorial distances in the schema and stores the result as an unordered map linking pairs of vertices \f$(i,j)\f$ where \f$i < j\f$ and the distance \f$d\f$ between them. 
    virtual void compute_distances(pair_index&) const;
    /// This method verifies if every vertex has at least one edge and so a positive valence, returning true in this case.
    bool positive_valence() const;
    /// This method computes a spanning tree of this 1-skeleton assumed to be connected, i.e. a minimal subset of edges of the 1-skeleton which keeps it connected, returning the number of edges used; the method's argument stores the edges as pairs of vertices.   
    int spanning_tree(std::vector<int>&) const;
    /// This method analyzes the number of connected components in the 1-skeleton, returning the number of distinct components and to which component each vertex belongs, the method's argument.
    int component_analysis(std::vector<int>&) const;
    /// This method returns true if there are no self-loops, every neighbour exists and these neighbour lists are mutually consistent.
    virtual bool consistent() const;
    /// This method returns the value of Schema::nvertex.
    int order() const;
    /// This method writes the value of the neighbour set to the second argument for the vertex whose index is the first argument.
    void get_neighbours(int,std::set<int>&) const;
    /// This method returns the degree of vertex in its unique argument.
    int degree(int) const;
  };

  inline int Schema::order() const 
  {
    return nvertex;
  }

  inline void Schema::get_neighbours(int n,std::set<int>& vx) const
  {
    vx = neighbours[n];
  }

  inline int Schema::add_vertex()
  {
    std::set<int> empty;
    int output = nvertex;

    neighbours.push_back(empty);
    nvertex++;

    return output;
  }

  inline bool Schema::positive_valence() const
  {
    // This method just checks if there are any vertices with no connections,
    // in which case it returns false, true otherwise
    for(int i=0; i<nvertex; ++i) {
      if (neighbours[i].empty()) return false;
    }
    return true;
  }

  inline int Schema::degree(int n) const
  {
    if (n < 0 || n >= nvertex) throw std::invalid_argument("The vertex argument in Schema::degree does not exist!");
    return (signed) neighbours[n].size();
  }

  inline bool Schema::connected(int n,int m) const
  {
    // A method to check if the vertices n and m share a direct 
    // connection
    if (n < 0 || n >= nvertex) throw std::invalid_argument("A vertex argument in Schema::connected does not exist!");
    if (m < 0 || m >= nvertex) throw std::invalid_argument("A vertex argument in Schema::connected does not exist!");
    if (n == m) throw std::invalid_argument("Vertex arguments are identical in the Schema::connected method!");

    // This edge exists...
    if (neighbours[n].count(m) > 0) return true;

    return false;
  }
}
#endif
