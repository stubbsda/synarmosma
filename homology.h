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
  std::string write() const;
  inline void set_method(METHOD m) {method = m;};
  inline void set_field(FIELD f) {field = f;};
  inline METHOD get_method() const {return method;};
  inline FIELD get_field() const {return field;};
  void initialize(FIELD,METHOD);
  void clear();
  void compute(const Nexus*);
  void append(const Group&);
  void serialize(std::ofstream&) const;
  void deserialize(std::ifstream&);
  void betti_numbers(std::vector<int>&) const;
};
#endif 
