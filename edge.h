#include "global.h"

#ifndef _edgeh
#define _edgeh

namespace SYNARMOSMA {
  /// A class representing an oriented 1-simplex or edge that connects two nodes in a graph or higher-dimensional simplicial complex.
  class Edge {
   protected:
    /// This integer property stores the lesser index of the two vertices connected 
    /// by this edge.
    int low = -1;
    /// This integer property stores the greater index of the two vertices connected 
    /// by this edge.
    int high = -1;
    /// This Boolean property is true if this edge is part of an n-cycle (n > 2) and 
    /// false otherwise.
    bool cyclic = false;
    /// This floating point property stores the edge's length and is this class's only 
    /// concession to geometry. If the edge's direction is considered causal in nature, 
    /// we can also assume that the length > 0 only when the direction is "disparate".
    double length = 0.0;
    /// This floating point property represents the flow of some quantity across the 
    /// edge and as such is positive for directed graphs (the flow must follow the edge's 
    /// innate orientation) but can be either sign in undirected graphs (positive then 
    /// means a flow from low to high, negative a flow in the opposite direction).
    double flow = 0.0;
    /// The capacity is always positive or zero and represents the maximum amount of 
    /// "mass" that can flow across the edge per unit time; in general we can imagine 
    /// that the capacity is in part a function of the length. 
    double capacity = 0.0;
    /// This property contains the edge's orientation: before (low => high), after (high => low) 
    /// or disparate (undirected).
    Relation direction = Relation::disparate;

    /// This method restores all the properties of the Edge instance to their default value. 
    void clear();
    /// This method returns the edge direction, considered from the ordering of the two arguments.
    inline Relation get_direction(int,int) const;
   public:
    /// The standard constructor that initializes all properties to its default value and exits.
    Edge();
    /// This constructor accepts arguments for the two vertices the edge connects and optionally its length and direction.
    Edge(int,int,double = 0.0,Relation = Relation::disparate);
    /// The standard copy constructor.
    Edge(const Edge&);
    /// The standard destructor.
    virtual ~Edge();
    /// The assignment operator, functions like the copy constructor. 
    Edge& operator =(const Edge&);
    /// This method writes the edge's properties to a binary disk file and returns the number of bytes written to the file.
    virtual int serialize(std::ofstream&) const;
    /// This method calls the clear() method on the instance and then reads the properties from a binary disk file and returns the number of bytes read.
    virtual int deserialize(std::ifstream&);
    /// This method copies the two vertices into the integer array argument, which is assumed to be of the form "int vx[2]".
    inline void get_vertices(int*) const;
    /// This method takes two arguments which are assumed to be the new values for the edge's vertices.
    inline void set_vertices(int,int);
    /// This method inverts the edge's orientation - if the edge is undirected it does nothing and returns false, otherwise it reverses the edge's direction and returns true.
    inline bool invert();
    /// This method returns 0 if the edge is undirected (disparate), +1 if low => high and -1 otherwise.
    inline int get_parity() const;
    friend class Graph;
    friend class Directed_Graph;
    friend class Logic_Graph;
  };

  void Edge::set_vertices(int u,int v) 
  {
    if (u == v) throw std::invalid_argument("The Edge class cannot be used for self-loops!");

    if (u < v) {
      low = u; high = v;
    }
    else {
      low = v; high = u;
    }
  }

  void Edge::get_vertices(int* vx) const
  {
    vx[0] = low; vx[1] = high;
  }

  int Edge::get_parity() const
  {
    if (direction == Relation::disparate) return 0;
    int output = (direction == Relation::before) ? 1 : -1;
    return output;
  }

  bool Edge::invert()
  {
    if (direction == Relation::disparate) return false;
    if (direction == Relation::before) {
      direction = Relation::after;
    }
    else {
      direction = Relation::before;
    }
    return true;
  }

  Relation Edge::get_direction(int u,int v) const
  {
    if ((u != low) && (u != high)) throw std::invalid_argument("The first argument of Edge::get_direction does not belong to the edge!"); 
    if ((v != low) && (v != high)) throw std::invalid_argument("The second argument of Edge::get_direction does not belong to the edge!"); 

    if (direction == Relation::disparate) return direction;
    if (u < v) {
      return direction;
    }
    else {
      if (direction == Relation::before) {
        return Relation::after;
      }
      else {
        return Relation::before;
      }
    }
  }
}
#endif

