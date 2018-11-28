#include "nexus.h"
#include "integer_matrix.h"

#ifndef _homologyh
#define _homologyh

namespace SYNARMOSMA {
  class Homology {
   public:
    enum class Field
    {
        int32,
        multiprecision,
        mod2
    };

    enum class Method
    {
        gap,
        native
    };

   private:
    std::vector<unsigned int> betti_number;
    std::vector<std::vector<unsigned int> > torsion;
    Method method = Method::native;
    Field field = Field::int32;

    void compute_integral_native(const Nexus*);
    void compute_native(const Nexus*);
    void compute_gap(const Nexus*);
  
   public:
    Homology();
    Homology(Field,Method);
    Homology(const Homology&);
    Homology& operator =(const Homology&);
    ~Homology();
    std::string write() const;
    inline void set_method(Method m) {method = m;};
    inline void set_field(Field f) {field = f;};
    inline Method get_method() const {return method;};
    inline Field get_field() const {return field;};
    inline void get_betti_numbers(std::vector<unsigned int>& output) const {output = betti_number;}; 
    void initialize(Field,Method);
    void clear();
    void compute(const Nexus*);
    int serialize(std::ofstream&) const;
    int deserialize(std::ifstream&);
  };
}
#endif 
