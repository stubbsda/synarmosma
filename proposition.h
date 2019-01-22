#include "random.h"

#ifndef _propositionh
#define _propositionh

namespace SYNARMOSMA {
  /// A class representing a logical proposition in conjunctive normal form, with a fixed number of atomic propositions per clause.
  class Proposition {
   protected:
    /// This vector of signed integers contains the proposition, stored according to the 
    /// following system: the vector contains an integral number of clauses and each clause  
    /// is of length 2*NP. It contains a non-negative integer (the index of the atomic 
    /// proposition) followed by 0 or 1, which indicates if the atomic proposition has a 
    /// NOT symbol. If the clause contains less than NP atomic propositions, then all of 
    /// its subsequent elements are set to -1. 
    std::vector<int> clause;
    /// This static property is an unsigned integer containing the maximum number of atomic 
    /// propositions in a clause.
    static const unsigned int NP = 5;

    /// This method accepts a value for the number of clauses to be created and the set of atomic propositions to be used and constructs the clause. 
    void initialize(unsigned int,const std::set<int>&);
    /// This method evaluates the truth value of the proposition, with the first argument a hash map of atomic propositions and their Boolean value and the second argument on output contains the indices of all the clauses that were false with this set of values for the atomic propositions.
    bool evaluate(const std::unordered_map<int,bool>&,std::set<unsigned int>&) const;
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
    /// This method evaluates the truth value of the proposition, given a vector of pairs that maps each atomic proposition to a Boolean value. 
    bool evaluate(const std::unordered_map<int,bool>&) const;
    /// This method uses Schöning's random restart algorithm to do a probabilistic test of the satisfiability of the proposition, i.e. does there exist a mapping of the atomic propositions such that this proposition is true.
    bool satisfiable() const;
    /// This method collects the distinct atoms contained in this proposition and writes them to the method's unique argument, returning the number of such atoms.
    unsigned int get_atoms(std::set<int>&) const;
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
    /// This method sets the Proposition::clause vector to be equal to the method's argument, after checking that its length is an even multiple of Proposition::NP.
    inline void set_clause(const std::vector<int>&);
    /// This method sets the argument to the Proposition::clause vector.
    inline void get_clause(std::vector<int>& c) const {c = clause;};
    /// This method returns the number of clauses in this proposition. 
    inline unsigned int get_size() const {return clause.size()/(2*NP);}; 
    /// This method returns the maximum number of atomic propositions in a clause.
    inline static unsigned int get_clause_size() {return NP;};
    /// This method overloads the ostream operator to do a pretty print of the proposition, with the atomic propositions written as "p[n]" where n >= 1 is an integer and using the standard logical connectives ∧, ∨ and ¬.
    friend std::ostream& operator <<(std::ostream&,const Proposition&);
    /// This method overloads the & operator to carry out a conjunction of the two propositions, thereby creating an output proposition whose size is the sum of the size of the two arguments.
    friend Proposition operator &(const Proposition&,const Proposition&);
  };

  void Proposition::set_clause(const std::vector<int>& c) {
    if (c.size()%(2*NP) != 0) throw std::invalid_argument("The argument of the Proposition::set_clause method has an illegal length!"); 
    clause = c;
  }
}
#endif
