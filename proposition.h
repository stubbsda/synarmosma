#include "global.h"

#ifndef _propositionh
#define _propositionh

namespace SYNARMOSMA {
  class Proposition {
   // A class to represent propositions in conjunctive normal form
   protected:
    std::vector<int> clause;

    // The number of atomic propositions in a clause
    static const int NP = 5;
   public:
    Proposition();
    Proposition(int);
    Proposition(const std::set<int>&);
    Proposition(int,const std::set<int>&);
    Proposition(const Proposition&);
    Proposition& operator =(const Proposition&);
    ~Proposition();
    void initialize(int,const std::set<int>&);
    bool evaluate(const std::vector<int>&,std::vector<int>&) const;
    bool evaluate(const std::vector<std::pair<int,bool> >&) const;
    bool satisfiable() const;
    void get_atoms(std::set<int>&) const;
    void set_atoms(const std::set<int>&);
    void mutate();
    void mutate(const std::set<int>&);
    inline void clear() {clause.clear();};
    void serialize(std::ofstream&) const;
    void deserialize(std::ifstream&); 
    inline static int get_clause_size() {return NP;};
    friend std::ostream& operator <<(std::ostream&,const Proposition&);
    friend Proposition operator &(const Proposition&,const Proposition&);
    friend class Logic_Graph;  
    friend class Propositional_System;
  };
}
#endif
