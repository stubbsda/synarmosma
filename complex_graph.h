#include "global.h"

#ifndef _cgraphh
#define _cgraphh

namespace SYNARMOSMA {
  class Complex_Graph {
   private:
    int n;
    std::vector<int>* neighbours;

    void compute_bridges() const;
   public:
    Complex_Graph();
    Complex_Graph(int);
    ~Complex_Graph();
    Complex_Graph(const Complex_Graph&);
    Complex_Graph& operator =(const Complex_Graph&);
    int get_candidates(std::vector<int>&) const;
    void add_edge(int,int);
    void contract(int,int,Complex_Graph&) const;
    void remove(int,int,Complex_Graph&) const;
  };
}
#endif

