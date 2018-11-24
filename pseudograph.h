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

    int DFS_bridge(int,int,int,int*,int*,hash_map&) const;
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
    int get_candidates(std::vector<int>&) const;
    /// This method calculates and returns the number of self-loops in this pseudograph.
    inline int get_loops() const;
    /// This method adds the edge joining the vertices indicated by the two arguments after checking that they lie within the appropriate range.
    void add_edge(int,int);
    void contract(int,int,Pseudograph*) const;
    void remove(int,int,Pseudograph*) const;
    friend class Graph;
  };

  int Pseudograph::get_loops() const
  {
    int i,sum = 0;
    for(i=0; i<nvertex; ++i) {
     sum += std::count(neighbours[i].begin(),neighbours[i].end(),i); 
    }
    return sum;
  }
}
#endif

