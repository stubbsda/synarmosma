#include "proposition.h"
#include "binary_matrix.h"
#include "directed_graph.h"

#ifndef _psystemh
#define _psystemh

namespace SYNARMOSMA {
  class Propositional_System {
   protected:
    std::vector<Proposition> theorems;
    unsigned int natom = 10;
    unsigned int nuniverse = 0;
    std::vector<boost::dynamic_bitset<> > truth;

    void compute_internal_logic();
    void initialize(int);
   public:
    Propositional_System();
    Propositional_System(int);
    Propositional_System(int,int);
    Propositional_System(int,const std::string&);
    Propositional_System(const Propositional_System&);
    Propositional_System& operator =(const Propositional_System&);
    ~Propositional_System();
    void clear();
    int serialize(std::ofstream&) const;
    int deserialize(std::ifstream&);
    inline int bit_count(int n) const {return truth[n].count();};
    int consistency(int,int,std::string&) const;
    bool implication(int,const std::vector<unsigned int>&) const;
    void compute_implication_graph(Directed_Graph*) const;
    friend class Logic_Graph;
  };
}
#endif 
