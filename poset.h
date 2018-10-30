#include "global.h"

#ifndef _poseth
#define _poseth

namespace SYNARMOSMA {
  /// A class representing a partially ordered set (poset), considered as non-negative integers, along with an ordering between certain pairs of elements. 
  class Poset {
   protected:
    /// This integer is simply the number of distinct elements in this 
    /// poset, which is assumed empty by default.
    int N = 0;
    /// This unordered_map from the Boost library contains the ordering on pairs of elements 
    /// of the poset; if x ~ y, then std::pair<x,y> is true.  
    boost::unordered_map<std::pair<int,int>,bool> order;

    void compute_width(int,int,std::set<int>&) const;
    int build_chain(std::vector<int>&,int) const;
   public:
    /// The default constructor which does nothing.
    Poset();
    /// This constructor sets the value of the property N to the argument but otherwise does nothing, i.e. no ordering on the elements.
    Poset(int);
    /// The copy constructor.
    Poset(const Poset&);
    /// The assignment operator.
    Poset& operator =(const Poset&);
    /// The destructor which does nothing.
    virtual ~Poset();
    /// This method sets N to 0 and empties the contents of the property order.
    virtual void clear();
    /// This method verifies that the poset's order property satisfies the axioms of an ordering, namely that is reflexive, anti-symmetric and transitive.
    virtual bool consistent() const;
    /// This method simply increments the property N by one. 
    inline void add_element() {N += 1;};
    /// This method tests if the element given by the argument is a sink, i.e. its posteriority is the empty set.
    bool sink(int) const;
    /// This method tests if the element given by the argument is a sink, i.e. its anteriority is the empty set.
    bool source(int) const;
    /// This method writes the instance properties to a binary disk file and returns the number of bytes written to the file.
    virtual int serialize(std::ofstream&) const;
    /// This method calls the clear() method on the instance and then reads the properties from a binary disk file and returns the number of bytes read.
    virtual int deserialize(std::ifstream&);
    /// This method computes the anteriority of an element x (the first argument), i.e. the set of all elements y of the poset such that y ~ x.
    void compute_anteriority(int,std::set<int>&) const;
    /// This method computes the posteriority of an element x (the first argument), i.e. the set of all elements y of the poset such that x ~ y.
    void compute_posteriority(int,std::set<int>&) const;
    /// This method accepts two elements x and y and adds the relation x ~ y to the order property of this poset, ensuring that the order continues to satisfy to be transitive; the method returns false if x = y or the relation already exists.
    bool set_order(int,int);
    /// This  method accepts two elements x and y and removes the relation x ~ y or y ~ x from the poset while respecting transtivity; the method returns false if x = y or the relation is disparate. 
    bool unset_order(int,int); 
    /// This method accepts two elements x and y and inverts their order, i.e. if x ~ y then y ~ x or vice-versa; the method returns false if x = y or the relation is disparate. 
    bool invert_order(int,int);
    /// This method constructs a random partial order on the set with a totality greater than or equal to the method's argument.
    void construct_ordering(double);
    void power_set(int);
    /// This method calculates the percentage of the N*(N-1)/2 distinct pairs (x,y) with x < y such that x ~ y, i.e. how close this partial order comes to being a total order.
    double totality() const;
    bool covered(int,int) const;
    int chain_number(int) const;
    void write_incastrature(const std::string&) const;
    /// This method accepts two elements x and y of the poset and returns "before" if x ~ y, "after" if x ~ y and "disparate" otherwise. 
    Relation get_order(int,int) const;
    friend class Directed_Graph;
  };
}
#endif

