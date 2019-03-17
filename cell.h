#include "global.h"

#ifndef _cellh
#define _cellh

namespace SYNARMOSMA {
  /// A class representing an n-simplex (n > 0) and which is used in the Schema class that represents an abstract simplicial complex.
  class Cell {
   protected:
    /// This property is a set containing all the 1+n vertices (by index) that this n-simplex is 
    /// built out of.  
    std::set<int> vertices;
    /// This integer set is a list of the index (in the Nexus class) of this n-simplex's n+1 faces that this n-simplex 
    /// posseses.
    std::set<int> entourage;
    /// This vector lists the vertex set of each of the n-simplex's 1+n faces, obtained by successively eliminating 
    /// one of the elements from Cell::vertices. 
    std::vector<std::set<int> > faces;

    /// This method computes the 1+n faces of this n-simplex and stores their vertex set in the vector Cell::faces.
    void calculate_faces();
   public:
    /// The default constructor, which does nothing.
    Cell();
    /// The standard copy constructor, which copies over the properties from the source instance of this class.
    Cell(const Cell&);
    /// This constructor builds an n-simplex from the 1+n vertices (0,1,...,n), where n is the constructor's positive argument.
    Cell(int);
    /// This is a specialized constructor to handle the common case of 1-simplices (edges) and accepts two arguments as the vertices connected by this edge.
    Cell(int,int);
    /// This is the normal constructor to use with this class, in which Cell::vertices is set to the argument and calculate_faces() is then called.
    Cell(const std::set<int>&);
    /// The destructor which in this case does nothing.
    virtual ~Cell();
    /// The standard overloaded assignment operator, which copies over the properties from the source instance of this class.
    Cell& operator =(const Cell&);
    /// This method checks if Cell::vertices contains the second argument; if so, it is swapped for the first argument, the calculate_faces() method is called and the method returns true, otherwise it return false.
    bool exchange(int,int);
    /// This method calls clear() and then creates a 1-simplex that connects the two vertices that are the method's arguments.
    virtual void initialize(int,int);
    /// This method calls clear() and then creates an n-simplex whose Cell::vertices property is the method's argument.
    virtual void initialize(const std::set<int>&);
    /// This method writes the instance properties to a binary disk file and returns the number of bytes written to the file.
    virtual int serialize(std::ofstream&) const;
    /// This method calls the clear() method on the instance and then reads the properties from a binary disk file and returns the number of bytes read.
    virtual int deserialize(std::ifstream&);
    /// This method clears the content of the Cell::vertices, Cell::entourage and Cell::faces properties.
    virtual void clear();
    /// This method checks if the argument is the vertex set of one of the faces of this n-simplex and returns true if this is the case.
    inline bool face(const std::set<int>&) const;
    /// This method writes the content of the Cell::vertices set to the method's argument, which is assumed to be long enough to store all 1+n vertices. 
    inline void get_vertices(int*) const;
    /// This method writes the content of the Cell::vertices set to the method's argument, as a vector.
    inline void get_vertices(std::vector<int>&) const;
    /// This method writes the content of the Cell::vertices set to the method's argument.
    inline void get_vertices(std::set<int>& v) const {v = vertices;};
    /// This method writes the content of the Cell::entourage set to the method's argument
    inline void get_entourage(std::set<int>& v) const {v = entourage;};
    /// This method writes the content of the Cell::faces vector to the method's argument.
    inline void get_faces(std::vector<std::set<int> >& v) const {v = faces;};
    /// This method returns the cardinality of Cell::vertices less one, i.e. the dimension of this simplex.
    inline int dimension() const {return (vertices.size() - 1);};
    /// This method returns true if the Cell::vertices contains the method's argument and false otherwise.
    inline bool contains(int) const;
    /// This method returns true if Cell::vertices is empty, false otherwise.
    inline bool empty() const {return vertices.empty();};
    /// This method tests if the dimension of its two arguments is the same (if not it returns zero), then computes the intersection of the Cell::vertices property of its two arguments and returns the size of the resulting set.
    friend inline int affinity(const Cell&,const Cell&);
    /// This operator returns true if the Cell::vertices property of both its arguments is the same, false otherwise.
    friend inline bool operator ==(const Cell&,const Cell&);
    /// This operator returns true if the Cell::vertices property of both its arguments isn't the same and false otherwise.
    friend inline bool operator !=(const Cell&,const Cell&);
    /// This operator returns true if the Cell::vertices property of the first argument is a subset of the Cell::vertices property of its second argument, false otherwise.
    friend inline bool operator <=(const Cell&,const Cell&);
    /// This operator returns true if the Cell::vertices property of the first argument is a proper subset of the Cell::vertices property of its second argument, false otherwise.
    friend inline bool operator <(const Cell&,const Cell&);
    friend class Nexus;
  };

  inline void Cell::get_vertices(std::vector<int>& v) const
  {
    std::set<int>::const_iterator it;

    v.clear();
    for(it=vertices.begin(); it!=vertices.end(); ++it) {
      v.push_back(*it);
    }
  }

  inline void Cell::get_vertices(int* v) const 
  {
    int n = 0;
    std::set<int>::const_iterator it;

    for(it=vertices.begin(); it!=vertices.end(); ++it) {
      v[n] = *it; n++;
    }
  }  

  inline bool Cell::contains(int v) const
  {
    bool output = false;
    if (vertices.count(v) > 0) output = true;
    return output;
  }

  bool Cell::face(const std::set<int>& f) const
  {
    if ((signed) f.size() != (dimension() - 1)) return false;
    for(int i=0; i<=dimension(); ++i) {
      if (f == faces[i]) return true;
    }
    return false;
  }

  inline int affinity(const Cell& c1,const Cell& c2)
  {
    int d = c1.dimension();

    if (d != c2.dimension()) return 0;

    std::set<int>::const_iterator it,jt;
    int i,j,nc = 0;

    for(it=c1.vertices.begin(); it!=c1.vertices.end(); ++it) {
      i = *it;
      for(jt=c2.vertices.begin(); jt!=c2.vertices.end(); ++jt) {
        j = *jt;
        if (i == j) nc++;
      }
    }

    return nc;
  }

  inline bool operator ==(const Cell& c1,const Cell& c2)
  {
    if (c1.vertices == c2.vertices) return true;
    return false;
  }

  inline bool operator !=(const Cell& c1,const Cell& c2)
  {
    if (c1 == c2) return false;
    return true;
  }

  inline bool operator <=(const Cell& c1,const Cell& c2)
  {
    if (c1.dimension() > c2.dimension()) return false;
    if (c1 == c2) return true;
    bool output = (c1 < c2) ? true : false;
    return output;
  }

  inline bool operator <(const Cell& c1,const Cell& c2)
  {
    if (c1.dimension() >= c2.dimension()) return false;
    // So c2 is bigger than c1, let's see if c1 is contained in c2
    return std::includes(c2.vertices.begin(),c2.vertices.end(),c1.vertices.begin(),c1.vertices.end());
  }
}
#endif

