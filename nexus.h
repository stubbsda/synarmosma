#include "cell.h"
#include "group.h"
#include "graph.h"
#include "matrix.h"

#ifndef _nexush
#define _nexush

// How to handle the base ZZ_p ?
enum FIELD {INT,ZZ,GF2};

enum METHOD {GAP,NATIVE};

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
  void compute_integral_homology(std::vector<Group>&,FIELD) const;
  void compute_homology_native(std::vector<Group>&,FIELD) const;
  void compute_homology(std::vector<Group>&,FIELD) const;
  void compute_homotopy(Group*) const;
  inline int get_dimension() const {return dimension;};
};
#endif


