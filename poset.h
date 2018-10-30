#include "global.h"

#ifndef _poseth
#define _poseth

namespace SYNARMOSMA {
  /// A class representing a partially ordered set (poset), considered as non-negative integers, along with a hash map that specifies the order (if any) between elements.
  class Poset {
   protected:
    int N = 0;
    boost::unordered_map<std::pair<int,int>,bool> order;

    void compute_width(int,int,std::set<int>&) const;
    int build_chain(std::vector<int>&,int) const;
   public:
    Poset();
    Poset(int);
    Poset(const Poset&);
    Poset& operator =(const Poset&);
    virtual ~Poset();
    virtual void clear();
    virtual bool consistent() const;
    inline void add_element() {N += 1;};
    bool sink(int) const;
    bool source(int) const;
    virtual int serialize(std::ofstream&) const;
    virtual int deserialize(std::ifstream&);
    void compute_anteriority(int,std::set<int>&) const;
    void compute_posteriority(int,std::set<int>&) const;
    bool set_order(int,int); 
    bool unset_order(int,int); 
    bool invert_order(int,int);
    void construct_ordering(double);
    void power_set(int);
    double totality() const;
    bool covered(int,int) const;
    int chain_number(int) const;
    void write_incastrature(const std::string&) const;
    Relation get_order(int,int) const;
    friend class Directed_Graph;
  };
}
#endif

