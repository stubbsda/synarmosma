#include "random.h"

#ifndef _propositionh
#define _propositionh

namespace SYNARMOSMA {
  /// A class representing a logical proposition in conjunctive normal form, with a fixed number of atomic propositions per clause.
  class Proposition {
   protected:
    /// This vector of signed integers contains the proposition, stored according to the 
    /// following system: the vector contains an integral number of clauses and each clause  
    /// is of length 2*Proposition::NP. It contains a non-negative integer (the index of 
    /// the atomic proposition) followed by 0 or 1, which indicates if the atomic proposition 
    /// has a NOT symbol. If the clause contains less than NP atomic propositions, then all 
    /// of its subsequent elements are set to -1. 
    std::vector<int> clause;

    /// This static property is an unsigned integer containing the maximum number of atomic 
    /// propositions in a clause.
    static const unsigned int NP = 5;
   public:
    /// The default constructor which does nothing.
    Proposition();
    Proposition(const std::set<int>&);
    Proposition(int,const std::set<int>&);
    Proposition(const Proposition&);
    Proposition& operator =(const Proposition&);
    ~Proposition();
    void initialize(unsigned int,const std::set<int>&);
    bool evaluate(const std::vector<int>&,std::vector<int>&) const;
    bool evaluate(const std::vector<std::pair<int,bool> >&) const;
    bool satisfiable() const;
    void get_atoms(std::set<int>&) const;
    void set_atoms(const std::set<int>&);
    void mutate();
    void mutate(const std::set<int>&);
    /// This method clears the Proposition::clause vector.
    inline void clear() {clause.clear();};
    int serialize(std::ofstream&) const;
    int deserialize(std::ifstream&);
    inline void set_clause(const std::vector<int>&);
    inline void get_clause(std::vector<int>& c) const {c = clause;};
    /// This method returns the number of clauses in this proposition. 
    inline unsigned int get_size() const {return clause.size()/(2*NP);}; 
    inline static unsigned int get_clause_size() {return NP;};
    friend std::ostream& operator <<(std::ostream&,const Proposition&);
    friend Proposition operator &(const Proposition&,const Proposition&);
    friend class Logic_Graph;  
    friend class Propositional_System;
  };

  void Proposition::set_clause(const std::vector<int>& c) {
    if (c.size()%(2*NP) != 0) throw std::invalid_argument("The argument of the Proposition::set_clause method has an illegal length!"); 
    clause = c;
  }
}
#endif
