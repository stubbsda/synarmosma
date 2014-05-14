#include "matrix.h"
#include "binary_matrix.h"
#include "group.h"
#include "nexus.h"

#ifndef _homologyh
#define _homologyh

// How to handle the base ZZ_p ?
enum FIELD {INT,ZZ,GF2};

enum METHOD {GAP,NATIVE};

class Homology {
 private:
  std::vector<Group> sequence;
  METHOD method;
  FIELD field;

  void compute_integral_native(const Nexus*);
  void compute_native(const Nexus*);
  void compute_gap(const Nexus*);
  
 public:
  Homology();
  Homology(FIELD,METHOD);
  ~Homology();
  unsigned int normalize_operator(const std::vector<int>*,int,int,std::vector<int>&) const;
  std::string write() const;
  void initialize(FIELD,METHOD);
  void clear();
  void compute(const Nexus*);
  void append(const Group&);
  void serialize(std::ofstream&) const;
  void deserialize(std::ifstream&);
  void betti_numbers(std::vector<int>&) const;
};
#endif 
