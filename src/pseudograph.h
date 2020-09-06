#include "global.h"

#ifndef _cgraphh
#define _cgraphh

namespace SYNARMOSMA {
  /// A class representing a pseudograph, i.e. a graph which permits multi-edges and self-loops.
  class Pseudograph {
   private:
    /// The number of vertices in the pseudograph, which is 
    /// also the size of the array of vectors that stores the 
    /// neighbour lists.
    int nvertex = 0;
    /// This array of integer vectors contains the list of 
    /// neighbours for each vertex of the pseudograph; we use 
    /// vectors rather than sets to allow for multi-edges.
    std::vector<int>* neighbours;

    /// This method performs a depth-first search for bridges in the pseudograph - the first two arguments are the vertices, the third is the depth, the next two integer arrays store which vertices have already been seen and finally a hash map to index the bridges found; the return value is the number of bridges found.
    int DFS_bridge(int,int,int,int*,int*,hash_map&) const;
    /// This method is the one normally called to carry out a complete analysis of bridges in the pseudograph and uses the DFS_bridge() method after initializing the arrays and the hash map; the return value is the total number of bridges in the pseudograph.
    int compute_bridges(hash_map&) const;
    /// This method frees the memory of the array Pseudograph::neighbours if Pseudograph::nvertex is greater than zero and then sets Pseudograph::nvertex to zero.
    void clear();
   public:
    /// The default constructor which does nothing.
    Pseudograph();
    /// A constructor which accepts the order of the pseudograph as its unique argument, setting the value of Pseudograph::nvertex and allocating the memory for Pseudograph::neighbours. 
    Pseudograph(int);
    /// The destructor which frees the memory of the array Pseudograph::neighbours if Pseudograph::nvertex is greater than zero.
    ~Pseudograph();
    /// The copy constructor that deletes the Pseudograph::neighbours array if necessary, allocates it after copying Pseudograph::nvertex and then copies the contents of this array.
    Pseudograph(const Pseudograph&);
    /// The overloaded assignment operator that behaves exactly like the copy constructor for this class.
    Pseudograph& operator =(const Pseudograph&);
    /// This method writes the instance properties to a binary disk file and returns the number of bytes written to the file.
    int serialize(std::ofstream&) const;
    /// This method calls the clear() method on the instance and then reads the properties from a binary disk file and returns the number of bytes read.
    int deserialize(std::ifstream&);
    /// This method uses the compute_bridges() method to compute all the bridges in the pseudograph - the output value of the method - and so assemble a list of candidate edges for contraction, i.e. those which are not bridges.
    int get_candidates(std::vector<int>&) const;
    /// This method simply returns the value of Pseudograph::nvertex, the order of the pseudograph.
    int order() const;
    /// This method calculates and returns the number of self-loops in this pseudograph.
    int get_loops() const;
    /// This method adds the edge joining the vertices indicated by the two arguments after checking that they lie within the appropriate range.
    void add_edge(int,int);
    /// This method fuses together two vertices (the first two arguments), removing any edges between them; it functions by replacing the greater of the two vertices by the lesser. The final argument is the output pseudograph in which the vertices have been fused.
    void contract(int,int,Pseudograph*) const;
    /// This method removes any edges between the two vertices that are the first arguments of the method, writing the new output pseudograph as the final argument.
    void remove(int,int,Pseudograph*) const;
  };

  inline int Pseudograph::order() const 
  {
    return nvertex;
  }

  inline int Pseudograph::get_loops() const
  {
    int i,sum = 0;
    for(i=0; i<nvertex; ++i) {
     sum += std::count(neighbours[i].begin(),neighbours[i].end(),i); 
    }
    return sum;
  }
}
#endif

