#include "global.h"

#ifndef _schemah
#define _schemah

namespace SYNARMOSMA {
  class Schema {
   protected:
    int nvertex;
    std::vector<std::set<int> > neighbours;

   public:
    Schema();
    Schema(int);
    virtual ~Schema();
    bool connected() const;
    virtual void clear();
    int add_vertex();
    bool positive_valence() const;
    int spanning_tree(std::vector<int>&) const;
    int component_analysis(std::vector<int>&) const;
    bool consistent() const;
    void components(std::vector<int>&,std::vector<int>&) const;
    inline int get_order() const {return nvertex;};
  };
}
#endif
