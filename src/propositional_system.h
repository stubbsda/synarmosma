#include "proposition.h"
#include "directed_graph.h"
#include "lattice.h"

#ifndef _psystemh
#define _psystemh

namespace SYNARMOSMA {
  /// A class representing a finite collection of instances of the Proposition class, using a shared set of atomic propositions.
  class Propositional_System {
   protected:
    /// The number of distinct possible truth value assignments 
    /// to these Propositional_System::natom atomic propositions, 
    /// \f$2^N\f$, where \f$N\f$ is the cardinality of Propositional_System::atoms.
    UINT64 nuniverse = 0;
    /// The set of atomic propositions which are used in at least 
    /// one of the propositions.
    std::set<int> atoms;
    /// The vector of propositions, each one of them 
    /// an instance of the Proposition class.
    std::vector<Proposition> theorems;

    /// This method computes the property the truth value of each element of Propositional_System::theorems for each possible assignment of values to the atomic propositions in Propositional_System::atoms, storing the result in the method's argument.
    void compute_internal_logic(std::vector<boost::dynamic_bitset<> >&); 
    /// This method has as its unique argument the number of theorems to be created, which are built randomly by the Proposition class.
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
    /// This method clears the Propositional_System::theorems and Propositional_System::atoms properties and sets Propositional_System::nuniverse to zero.
    void clear();
    /// This method takes two elements of Propositional_System::theorems, the vector of dynamic bitsets computed in the compute_internal_logic() method and a Boolean operator specified by a string among "AND", "OR" and "XOR", and computes the number of logical universes in which these two propositions are true according to the Boolean operator. 
    unsigned int consistency(unsigned int,unsigned int,const std::string&,const std::vector<boost::dynamic_bitset<> >&) const;
    /// This method tests whether or not a given proposition is implied by a set of other propositions, in which case it returns true. The first argument is the element of Propositional_System::theorems being tested, the second is the set of axioms, which are other elements of this same vector and the third if the vector of dynamic bitsets calculated by the compute_internal_logic() method. A proposition q implies another p if (q or !p) is always true. 
    bool implication(unsigned int,const std::vector<unsigned int>&,const std::vector<boost::dynamic_bitset<> >&) const;
    /// This method constructs the directed graph of implications among the elements of Propositional_System::theorems; each proposition is a vertex and if one proposition p implies another q, then the edge p -> q is added to the graph. It requires the vector of dynamic bitsets calculated by the compute_internal_logic() method.
    void compute_implication_graph(Directed_Graph*,const std::vector<boost::dynamic_bitset<> >&) const;
    /// This method writes the instance properties to a binary disk file and returns the number of bytes written to the file.
    int serialize(std::ofstream&) const;
    /// This method calls the clear() method on the instance and then reads the properties from a binary disk file and returns the number of bytes read.
    int deserialize(std::ifstream&);
    /// This method appends a new element to the vector Propositional_System::theorems, then checks to see if any new atomic propositions have to be added, recomputing Propositional_System::nuniverse if so, and then returns the index of this new theorem.
    unsigned int add_theorem(const std::set<int>&);
    /// This method removes the proposition from the Propositional_System::theorems property, corresponding to the proposition whose index is the method's unique argument.
    void drop_theorem(unsigned int);
    /// This method returns the cardinality of the set Propositional_System::atoms.
    unsigned int get_number_atoms() const;
    /// This method sets the second argument to the set of atomic propositions for the proposition whose index is the first argument. 
    void get_atoms(unsigned int,std::set<int>&) const;
    /// This method returns the truth value of a given proposition in Propositional_System::theorems (whose index is the first argumemt), when its atoms are assigned the truth values specified in the second argument. 
    bool evaluate(unsigned int,const std::unordered_map<int,bool>&) const;
    /// This method computes the lattice corresponding to the subset relation among the set of atomic propositions for each element of the Propositional_System::theorems property, which is then written to the method's argument. The return value is that of the consistent method called on the argument after the lattice has been constructed.
    bool compute_propositional_lattice(Lattice*) const;
    friend class Logic_Graph;
  };

  inline unsigned int Propositional_System::get_number_atoms() const 
  {
    return atoms.size();
  }

  inline void Propositional_System::drop_theorem(unsigned int n)
  {
    if (n >= theorems.size()) throw std::invalid_argument("Illegal theorem index in the Propositional_System::drop_theorem method!");

    theorems.erase(theorems.begin() + n);
  }

  inline void Propositional_System::get_atoms(unsigned int n,std::set<int>& atoms) const 
  {
    if (n >= theorems.size()) throw std::invalid_argument("Illegal theorem index in the Propositional_System::get_atoms method!");
    theorems[n].get_atoms(atoms);
  }

  inline bool Propositional_System::evaluate(unsigned int n,const std::unordered_map<int,bool>& atoms) const
  {
    if (n >= theorems.size()) throw std::invalid_argument("Illegal theorem index in the Propositional_System::evaluate method!");
    return theorems[n].evaluate(atoms);
  }
}
#endif 
