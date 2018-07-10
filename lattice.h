#include "poset.h"

#ifndef _latticeh
#define _latticeh

namespace SYNARMOSMA {
  class Lattice: public Poset {
   protected:
    std::set<int> atoms;
    bool atomic;
    int null;
    int unity;
     
    void compute_atoms(); 
    void compute_bounds();
    void initialize();
   public:
    Lattice();
    Lattice(int);
    Lattice(const Lattice&);
    ~Lattice() override;
    Lattice& operator =(const Lattice&);
    int meet(int,int) const;
    int join(int,int) const;
    int serialize(std::ofstream&) const override;
    int deserialize(std::ifstream&) override;
    void clear() override;
    inline void add_element() {N += 1;};
    bool consistent() const override;
  };
}
#endif

