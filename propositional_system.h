#include "proposition.h"
#include <boost/dynamic_bitset.hpp>

#ifndef _psystemh
#define _psystemh

class Propositional_System {
 private:
  std::vector<Proposition> theorems;
  unsigned int natom;
  unsigned int nuniverse;
  bool** logical_universe;
  std::vector<boost::dynamic_bitset<> > truth;

  void compute_internal_logic();
  void build_logical_universe();
  void set_default_values();
  void initialize(unsigned int);
  
 public:
  Propositional_System();
  Propositional_System(unsigned int);
  Propositional_System(unsigned int,unsigned int);
  Propositional_System(unsigned int,const char*);
  ~Propositional_System();
  void read(const char*,unsigned int);
  void read(const char*);
  void write(const char*,unsigned int);
  void write(const char*);
  unsigned int bit_count(unsigned int) const;
  unsigned int consistency(unsigned int,unsigned int,const std::string&) const;
  bool implication(const Proposition&) const;
  bool implication(unsigned int,const std::vector<unsigned int>&) const;
  void compute_pairs(unsigned int*);
  friend class Logic_Graph;
};
#endif 
