#include "global.h"

#ifndef _poseth
#define _poseth

namespace SYNARMOSMA {
  class Poset {
   protected:
    unsigned int N;
    std::unordered_map<std::pair<unsigned int,unsigned int>,bool> order;

    void compute_width(unsigned int,unsigned int,std::set<unsigned int>&) const;
    unsigned int build_chain(std::vector<unsigned int>&,unsigned int) const;
    virtual void serialize(std::ofstream&) const;
    virtual void deserialize(std::ifstream&);
   public:
    Poset();
    Poset(unsigned int);
    Poset(const Poset&);
    Poset& operator =(const Poset&);
    virtual ~Poset();
    virtual void clear();
    virtual bool consistent() const;
    inline void add_element() {N += 1;};
    bool sink(unsigned int) const;
    bool source(unsigned int) const;
    void compute_anteriority(unsigned int,std::set<unsigned int>&) const;
    void compute_posteriority(unsigned int,std::set<unsigned int>&) const;
    bool set_order(unsigned int,unsigned int); 
    bool unset_order(unsigned int,unsigned int); 
    bool invert_order(unsigned int,unsigned int);
    void construct_ordering(double);
    void power_set(unsigned int);
    double totality() const;
    bool covered(unsigned int,unsigned int) const;
    unsigned int chain_number(unsigned int) const;
    void write_incastrature(const std::string&) const;
    RELATION get_order(unsigned int,unsigned int) const;
    friend class Directed_Graph;
  };
}
#endif

