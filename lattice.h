#include "global.h"

#ifndef _latticeh
#define _latticeh

namespace SYNARMOSMA {
  enum RELATION
  {
      BEFORE,
      AFTER,
      INCOMPARABLE
  };

  class Lattice {
   private:
    int N;
    boost::unordered_map<std::string,RELATION> order;
 
   public:
    Lattice();
    Lattice(int);
    Lattice(const Lattice&);
    ~Lattice();
    void clear();
    void add_vertex();
    RELATION get_relation(int,int) const;
    bool consistent() const;
  };
}
#endif

