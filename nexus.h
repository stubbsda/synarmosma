#include "cell.h"
#include "graph.h"

#ifndef _nexush
#define _nexush

namespace SYNARMOSMA {
  class Nexus : public Schema {
   protected:
    int dimension;
    std::vector<Cell>* elements;
    hash_map* index_table;

   public:
    Nexus();
    Nexus(int);
    Nexus(const Nexus&);
    Nexus& operator =(const Nexus&);
    virtual ~Nexus();
    virtual int serialize(std::ofstream&) const;
    virtual int deserialize(std::ifstream&);
    bool orientable() const;
    bool pseudomanifold(bool*) const;
    void surface_construction(int);
    virtual void clear();
    void initialize(int);
    void initialize(int,int);
    void paste(const std::set<int>&);
    void paste(const Cell&);
    void regularization();
    void ascend(int,int,std::vector<Cell>&) const;
    void star(const std::set<std::set<int> >&,std::vector<Cell>*) const;
    void link(const std::set<std::set<int> >&,std::vector<Cell>*) const;
    void closure(const std::set<std::set<int> >&,Nexus*,int*) const;
    void compute_entourages();
    void compute_neighbours();
    inline int get_dimension() const {return dimension;};
    friend class Homology;
    friend class Homotopy;
  };
}
#endif


