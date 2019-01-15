#include "cell.h"
#include "graph.h"

#ifndef _nexush
#define _nexush

namespace SYNARMOSMA {
  /// A class representing an abstract simplicial complex, using the Cell class to store higher-dimensional elements.
  class Nexus: public Schema {
   protected:
    int dimension = -1;
    std::vector<Cell>* elements;
    hash_map* index_table;

   public:
    Nexus();
    Nexus(int);
    Nexus(const Nexus&);
    Nexus& operator =(const Nexus&);
    ~Nexus() override;
    int serialize(std::ofstream&) const override;
    int deserialize(std::ifstream&) override;
    bool orientable() const;
    bool pseudomanifold(bool*) const;
    void surface_construction(std::string&);
    void clear() override;
    void initialize(int);
    void initialize(int,int);
    inline bool paste(const std::set<int>& vx) {return paste(Cell(vx));};
    bool paste(const Cell&);
    void regularization();
    void ascend(int,int,std::vector<Cell>&) const;
    void star(const std::set<std::set<int> >&,std::vector<Cell>*) const;
    void link(const std::set<std::set<int> >&,std::vector<Cell>*) const;
    void closure(const std::set<std::set<int> >&,Nexus*,int*) const;
    void compute_entourages();
    void compute_neighbours();
    inline int get_dimension() const {return dimension;};
    inline unsigned int get_length(int n) const {return elements[n].size();};
    inline int get_index(int D,std::set<int>& S) const {hash_map::const_iterator qt = index_table[D].find(S); return qt->second;};
    inline void get_elements(int D,int n,int* vx) const {elements[D][n].get_vertices(vx);};
    inline void get_elements(int D,int n,std::set<int>& vx) const {elements[D][n].get_vertices(vx);};
    inline void get_entourage(int D,int n,std::set<int>& vx) const {elements[D][n].get_entourage(vx);};
  };
}
#endif


