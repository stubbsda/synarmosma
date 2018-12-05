#include "proposition.h"
#include "directed_graph.h"

#ifndef _psystemh
#define _psystemh

namespace SYNARMOSMA {
  /// A class representing a finite collection of instances of the Proposition class, using a shared set of atomic propositions.
  class Propositional_System {
   protected:
    /// The number of distinct possible truth value assignments 
    /// to these Propositional_System::natom atomic propositions, 
    /// \f$2^N\f$, where \f$N\f$ is the cardinality of Propositional_System::atoms.
    unsigned int nuniverse = 0;
    /// The set of atomic propositions which are used in at least 
    /// one of the propositions.
    std::set<int> atoms;
    /// A vector containing all of the possible truth value 
    /// assignments for the atomic propositions, stored in each 
    /// case a bitset whose length is the cardinality of Propositional_System::atoms.
    std::vector<boost::dynamic_bitset<> > truth;
    /// The vector of propositions, each one of them 
    /// an instance of the Proposition class.
    std::vector<Proposition> theorems;

    /// This method computes the property Propositional_System::truth, first calculating the value of Propositional_System::nuniverse.
    void compute_internal_logic();
    /// This method has as its unique argument the number of theorems to be created, which are built randomly by the Proposition class, and then calls compute_internal_logic().
    void initialize(unsigned int);
   public:
    /// The default constructor which does nothing.
    Propositional_System();
    /// This is the standard constructor for this class, with two arguments equal to the number of atomic propositions and the number of theorems, and which calls initialize() after setting these values.
    Propositional_System(const std::set<int>&,unsigned int);
    /// The standard copy constructor that simply copies over the instance's four properties to the target instance.
    Propositional_System(const Propositional_System&);
    /// Overloaded assignement operator that copies over the instance's four properties to the target instance.
    Propositional_System& operator =(const Propositional_System&);
    /// The destructor which in the case of this class has nothing to do.
    ~Propositional_System();
    /// This method clears the Propositional_System::theorems and Propositional_System::truth properties and sets the other properties to zero.
    void clear();
    /// This method takes two elements of Propositional_System::theorems along with a Boolean operator among "AND", "OR" and "XOR", and computes the number of logical universes in which these two propositions are true according to the Boolean operator. 
    unsigned int consistency(unsigned int,unsigned int,std::string&) const;
    /// This method tests whether or not a given proposition is implied by a set of other propositions, in which case it returns true. The first argument is the element of Propositional_System::theorems being tested, the second is the set of axioms, which are other elements of this same vector. A proposition q implies another p if (q or !p) is always true. 
    bool implication(unsigned int,const std::vector<unsigned int>&) const;
    /// This method constructs the directed graph of implications among the elements of Propositional_System::theorems; each proposition is a vertex and if one proposition p implies another q, then the edge p -> q is added to the graph.
    void compute_implication_graph(Directed_Graph*) const;
    /// This method writes the instance properties to a binary disk file and returns the number of bytes written to the file.
    int serialize(std::ofstream&) const;
    /// This method calls the clear() method on the instance and then reads the properties from a binary disk file and returns the number of bytes read.
    int deserialize(std::ifstream&);
    /// This method appends a new element to the vector Propositional_System::theorems, then checks to see if any new atomic propositions have to be added (in which case it calls compute_internal_logic()) and then returns the index of this new theorem.
    unsigned int add_theorem(const std::set<int>&);
    /// This method removes the proposition from the Propositional_System::theorems and Propositional_System::truth properties, corresponding to the proposition whose index is the method's unique argument.
    bool drop_theorem(unsigned int);
    /// This method returns the cardinality of the set Propositional_System::atoms.
    inline unsigned int get_number_atoms() const {return atoms.size();};
    /// This method returns the number of logical universes in which a particular element of Propositional_System::theorems (specified by the method's unique argument) is true.
    inline unsigned int bit_count(unsigned int) const;
    /// This method sets the second argument to the set of atomic propositions for the proposition whose index is the first argument. 
    inline void get_atoms(unsigned int i,std::set<int>& atoms) const;
  };

  unsigned int Propositional_System::bit_count(unsigned int n) const
  {
    if (n >= theorems.size()) throw std::invalid_argument("Illegal theorem index in Propositional_System class!");
    return truth[n].count();
  }

  void Propositional_System::get_atoms(unsigned int n,std::set<int>& atoms) const 
  {
    if (n >= theorems.size()) throw std::invalid_argument("Illegal theorem index in Propositional_System class!");
    theorems[n].get_atoms(atoms);
  }
}
#endif 
