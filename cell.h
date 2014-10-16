#include "global.h"

#ifndef _cellh
#define _cellh

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
  void initialize(const std::set<int>&);
  void calculate_faces();
  virtual void clear();
  inline int dimension() const {return (vertices.size() - 1);};
  bool exchange(int,int);
  bool face(const std::set<int>&) const;
  double dimensional_stress(int) const;
  void get_vertices(int*) const;
  void get_vertices(std::vector<int>&) const;
  void get_vertices(std::set<int>&) const;
  void get_faces(std::vector<Cell>&) const;
  void serialize(std::ofstream&) const;
  virtual void deserialize(std::ifstream&);
  bool contains(int) const;
  inline bool empty() const {return vertices.empty();};
  friend int affinity(const Cell&,const Cell&);
  friend bool operator ==(const Cell&,const Cell&);
  friend bool operator !=(const Cell&,const Cell&);
  friend bool operator <=(const Cell&,const Cell&);
  friend bool operator <(const Cell&,const Cell&);
  friend Cell operator ^(const Cell&,const Cell&);
  friend std::ostream& operator<< (std::ostream&,const Cell&);
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
  std::set<int>::const_iterator it = std::find(vertices.begin(),vertices.end(),v);
  bool output = (it == vertices.end()) ? false : true;
  return output;
}

inline void Cell::get_vertices(std::set<int>& v) const
{
  v = vertices;
}

inline int affinity(const Cell& s1,const Cell& s2)
{
  int d = s1.dimension();

  if (d != s2.dimension()) return 0;

  std::set<int>::const_iterator it,jt;
  int i,j,nc = 0;

  for(it=s1.vertices.begin(); it!=s1.vertices.end(); ++it) {
    i = *it;
    for(jt=s2.vertices.begin(); jt!=s2.vertices.end(); ++jt) {
      j = *jt;
      if (i == j) nc++;
    }
  }

  return nc;
}
#endif

