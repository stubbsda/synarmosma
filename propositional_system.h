#include "proposition.h"
#include "binary_matrix.h"
#include "directed_graph.h"

#ifndef _psystemh
#define _psystemh

namespace SYNARMOSMA {
  /// A class representing a finite collection of instances of the Proposition class, using a shared set of atomic propositions.
  class Propositional_System {
   protected:
    /// The number of atomic propositions in the entire set of 
    /// propositions and which are presumed to be indexed from 
    /// 0 to Propositional_System::natom-1.
    unsigned int natom = 0;
    /// The number of distinct possible truth value assignments 
    /// to these Propositional_System::natom atomic propositions, 
    /// \f$2^N\f$, where \f$N\f$ is Propositional_System::natom.
    unsigned int nuniverse = 0;
    /// A vector containing all of the possible truth value 
    /// assignments for the atomic propositions, stored in each 
    /// case a bitset of length Propositional_System::natom.
    std::vector<boost::dynamic_bitset<> > truth;
    /// The vector of propositions, each one of them 
    /// an instance of the Proposition class.
    std::vector<Proposition> theorems;

    void compute_internal_logic();
    void initialize(unsigned int);
   public:
    /// The default constructor which does nothing.
    Propositional_System();
    Propositional_System(unsigned int,unsigned int);
    Propositional_System(const Propositional_System&);
    Propositional_System& operator =(const Propositional_System&);
    ~Propositional_System();
    /// This method clears the Propositional_System::theorems and Propositional_System::truth properties and sets the other properties to zero.
    void clear();
    unsigned int consistency(unsigned int,unsigned int,std::string&) const;
    bool implication(unsigned int,const std::vector<unsigned int>&) const;
    void compute_implication_graph(Directed_Graph*) const;
    /// This method writes the instance properties to a binary disk file and returns the number of bytes written to the file.
    int serialize(std::ofstream&) const;
    /// This method calls the clear() method on the instance and then reads the properties from a binary disk file and returns the number of bytes read.
    int deserialize(std::ifstream&);
    unsigned int add_theorem(const std::set<int>&);
    bool drop_theorem(unsigned int);
    inline unsigned int get_number_atoms() const {return natom;};
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
