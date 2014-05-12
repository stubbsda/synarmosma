#include "matrix.h"
#include "binary_matrix.h"
#include "group.h"

#ifndef _homologyh
#define _homologyh

enum FIELD {INTEGER,ZZ,GF2};

enum METHOD {GAP,NATIVE};

class Homology {
 private:
  std::vector<Group> sequence;
  METHOD method;
  FIELD field;
  
 public:
  Homology();
  Homology(FIELD,METHOD);
  ~Homology();
  unsigned int normalize_operator(const std::vector<int>*,int,int,std::vector<int>&) const;
  void clear();
  void append(const Group&);
  void betti_numbers(std::vector<int>&) const;
};
#endif 
