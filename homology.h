#include "nexus.h"

#ifndef _homologyh
#define _homologyh

namespace SYNARMOSMA {
  class Homology {
   public:
    enum FIELD
    {
        INT,
        ZZ,
        GF2
    };

    enum METHOD
    {
        GAP,
        NATIVE
    };

   private:
    std::vector<unsigned int> betti_number;
    std::vector<std::vector<unsigned int> > torsion;
    METHOD method;
    FIELD field;

    void compute_integral_native(const Nexus*);
    void compute_native(const Nexus*);
    void compute_gap(const Nexus*);
  
   public:
    Homology();
    Homology(FIELD,METHOD);
    Homology(const Homology&);
    Homology& operator =(const Homology&);
    ~Homology();
    std::string write() const;
    inline void set_method(METHOD m) {method = m;};
    inline void set_field(FIELD f) {field = f;};
    inline METHOD get_method() const {return method;};
    inline FIELD get_field() const {return field;};
    inline void get_betti_numbers(std::vector<unsigned int>& output) const {output = betti_number;}; 
    void initialize(FIELD,METHOD);
    void clear();
    void compute(const Nexus*);
    int serialize(std::ofstream&) const;
    int deserialize(std::ifstream&);
  };
}
#endif 
