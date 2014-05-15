#include "cell.h"
#include "graph.h"

#ifndef _nexush
#define _nexush

class Nexus : public Schema {
 private:
  int dimension;
  std::vector<Cell>* elements;
  hash_map* index_table;

 public:
  Nexus();
  Nexus(int);
  virtual ~Nexus();
  bool orientable() const;
  bool pseudomanifold(bool*) const;
  void assemble();
  virtual void clear();
  void initialize(int);
  void initialize(int,int);
  void paste(const std::set<int>&);
  void regularization();
  int size() const;
  void ascend(int,int,std::vector<Cell>&) const;
  void star(const std::set<std::string>&,std::vector<Cell>*) const;
  void link(const std::set<std::string>&,std::vector<Cell>*) const;
  void closure(const std::set<std::string>&,Nexus*,int*) const;
  void compute_entourages();
  void compute_neighbours();
  inline int get_dimension() const {return dimension;};
  friend class Homology;
  friend class Homotopy;
};
#endif


