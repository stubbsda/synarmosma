#include "global.h"

#ifndef _propositionh
#define _propositionh

class Proposition {
 private:
  std::vector<int> clause;

 public:
  Proposition();
  Proposition(int);
  Proposition(const std::set<int>&);
  Proposition(int,const std::set<int>&);
  Proposition(const Proposition&);
  Proposition& operator =(const Proposition&);
  Proposition& operator *(const Proposition&);
  ~Proposition();
  void initialize(int,const std::set<int>&);
  bool evaluate(const std::vector<int>&,std::vector<int>&) const;
  bool evaluate(const bool*) const;
  bool satisfiable() const;
  void atoms(std::set<int>&) const;
  inline void clear() {clause.clear();};
  void serialize(std::ofstream&) const;
  void deserialize(std::ifstream&); 
  friend std::ostream& operator <<(std::ostream&,const Proposition&);
  friend class Logic_Graph;  
  friend class Propositional_System;
};
#endif
