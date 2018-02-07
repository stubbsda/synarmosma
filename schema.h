#include "global.h"

#ifndef _schemah
#define _schemah

namespace SYNARMOSMA {
  class Schema {
   protected:
    int nvertex;
    std::vector<std::set<int> > neighbours;

   public:
    Schema();
    Schema(int);
    Schema(const Schema&);
    Schema& operator =(const Schema&);
    virtual ~Schema();
    virtual int serialize(std::ofstream&) const;
    virtual int deserialize(std::ifstream&);
    bool connected() const;
    inline bool connected(int,int) const;
    virtual void clear();
    inline int add_vertex();
    virtual bool add_edge(int,int,double = 0.0);
    virtual bool drop_edge(int,int);
    virtual int distance(int,int) const;
    virtual void compute_distances(edge_hash&) const;
    bool positive_valence() const;
    int spanning_tree(std::vector<int>&) const;
    int component_analysis(std::vector<int>&) const;
    virtual bool consistent() const;
    void components(std::vector<int>&,std::vector<int>&) const;
    inline int get_order() const {return nvertex;};
  };

  int Schema::add_vertex()
  {
    std::set<int> empty;
    neighbours.push_back(empty);
    nvertex++;
    return nvertex-1;
  }

  bool Schema::connected(int n,int m) const
  {
    // A method to check if the vertices n and m share a direct 
    // connection
#ifdef DEBUG
    assert(n >= 0 && n < nvertex);
    assert(m >= 0 && m < nvertex);
    if (n == m) return false;
#endif

    if (neighbours[n].count(m) == 0) {
      // This edge doesn't exist...
#ifdef DEBUG
      assert(neighbours[m].count(n) == 0);
#endif
      return false;
    }
    return true;
  }
}
#endif
