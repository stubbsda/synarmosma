#include "poset.h"

#ifndef _latticeh
#define _latticeh

namespace SYNARMOSMA {
  class Lattice: public Poset {
   protected:
    std::set<unsigned int> atoms;
    bool atomic;
    unsigned int null;
    unsigned int unity;
     
    void compute_atoms(); 
    void compute_bounds();
    void initialize();
   public:
    Lattice();
    Lattice(unsigned int);
    Lattice(const Lattice&);
    ~Lattice() override;
    Lattice& operator =(const Lattice&);
    unsigned int meet(unsigned int,unsigned int) const;
    unsigned int join(unsigned int,unsigned int) const;
    int serialize(std::ofstream&) const override;
    int deserialize(std::ifstream&) override;
    void clear() override;
    inline void add_element() {N += 1;};
    bool consistent() const override;
  };
}
#endif

