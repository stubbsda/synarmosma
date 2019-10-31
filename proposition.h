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

    void initialize(unsigned int,const std::set<int>&);
   public:
    /// The default constructor which does nothing.
    Proposition();
    /// This constructor accepts a set containing the atomic propositions which should be used in this proposition; the number of clauses is chosen to a random integer between 1 and the cardinality of the input divided by Proposition::NP.
    Proposition(const std::set<int>&);
    /// This constructor has two arguments, the number of clauses and a set of atomic propositions that should be used to build the proposition. 
    Proposition(unsigned int,const std::set<int>&);
    /// The standard copy constructor that copies over the value of the Proposition::clause property.
    Proposition(const Proposition&);
    /// The overloaded assignment operator that copies over the value of the Proposition::clause property.
    Proposition& operator =(const Proposition&);
    /// The destructor which does nothing for this class.
    ~Proposition();
    bool evaluate(const std::vector<int>&,std::vector<int>&) const;
    bool evaluate(const std::vector<std::pair<int,bool> >&) const;
    bool satisfiable() const;
    /// This method collects the distinct atoms contained in this proposition and writes them to the method's unique argument.
    void get_atoms(std::set<int>&) const;
    /// This method functions in exactly the same manner as the constructor with the same argument - a random number of clauses is generated and the atomic propositions are chosen from among the input set.
    void set_atoms(const std::set<int>&);
    /// This method retains the same number of clauses and the same set of atomic propositions but otherwise completely regenerates the proposition.
    void mutate();
    /// This method clears the Proposition::clause vector.
    inline void clear() {clause.clear();};
    /// This method writes the instance properties to a binary disk file and returns the number of bytes written to the file.
    int serialize(std::ofstream&) const;
    /// This method calls the clear() method on the instance and then reads the properties from a binary disk file and returns the number of bytes read.
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
