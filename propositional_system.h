#include "proposition.h"
#include "binary_matrix.h"
#include "directed_graph.h"

#ifndef _psystemh
#define _psystemh

namespace SYNARMOSMA {
  class Propositional_System {
   protected:
    std::vector<Proposition> theorems;
    unsigned int natom;
    unsigned int nuniverse;
    std::vector<boost::dynamic_bitset<> > truth;

    void compute_internal_logic();
    void set_default_values();
    void initialize(unsigned int);
  
   public:
    Propositional_System();
    Propositional_System(unsigned int);
    Propositional_System(unsigned int,unsigned int);
    Propositional_System(unsigned int,const char*);
    Propositional_System(const Propositional_System&);
    Propositional_System& operator =(const Propositional_System&);
    ~Propositional_System();
    void clear();
    int serialize(std::ofstream&) const;
    int deserialize(std::ifstream&);
    unsigned int bit_count(unsigned int) const;
    unsigned int consistency(unsigned int,unsigned int,const std::string&) const;
    bool implication(unsigned int,const std::vector<unsigned int>&) const;
    void compute_implication_graph(Directed_Graph*) const;
    friend class Logic_Graph;
  };
}
#endif 
