#include "global.h"

#ifndef _cellh
#define _cellh

namespace SYNARMOSMA {
  class Cell {
   protected:
    std::set<int> vertices;
    std::set<int> entourage;
    std::vector<std::set<int> > faces;

   public:
    Cell();
    Cell(const Cell&);
    Cell(int);
    Cell(int,int);
    Cell(const std::set<int>&);
    virtual ~Cell();
    Cell& operator =(const Cell&);
    void initialize(int,int);
    virtual void initialize(const std::set<int>&);
    void calculate_faces();
    virtual void clear();
    inline int dimension() const {return (vertices.size() - 1);};
    bool exchange(int,int);
    bool face(const std::set<int>&) const;
    void inline get_vertices(int*) const;
    void inline get_vertices(std::vector<int>&) const;
    void inline get_vertices(std::set<int>& v) const {v = vertices;};
    void get_faces(std::vector<Cell>&) const;
    virtual int serialize(std::ofstream&) const;
    virtual int deserialize(std::ifstream&);
    inline bool contains(int) const;
    inline bool empty() const {return vertices.empty();};
    friend inline int affinity(const Cell&,const Cell&);
    friend inline bool operator ==(const Cell&,const Cell&);
    friend inline bool operator !=(const Cell&,const Cell&);
    friend inline bool operator <=(const Cell&,const Cell&);
    friend inline bool operator <(const Cell&,const Cell&);
    friend class Nexus;
    friend class Homology;
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
    bool output = (vertices.count(v) == 0) ? false : true;
    return output;
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
    if (c2.dimension() < c1.dimension()) return false;
    if (c1 == c2) return true;
    bool output = (c1 < c2) ? true : false;
    return output;
  }

  inline bool operator <(const Cell& c1,const Cell& c2)
  {
    if (c2.dimension() <= c1.dimension()) return false;
    // So c2 is bigger than c1, let's see if c1 is contained
    // in c2
    return std::includes(c2.vertices.begin(),c2.vertices.end(),c1.vertices.begin(),c1.vertices.end());
  }
}
#endif

