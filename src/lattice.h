#include "poset.h"

#ifndef _latticeh
#define _latticeh

namespace SYNARMOSMA {
  /// A class representing a lattice, that is an instance of the Poset class in which every two elements have a unique supremum and a unique infimum.
  class Lattice: public Poset {
   protected:
    /// This property stores the atoms of the lattice, i.e. those elements \f$0\sim x\f$ such that there does not exist any \f$y\f$ satisfying \f$0\sim y\sim x\f$. 
    std::set<int> atoms;
    /// A lattice is atomic if for every non-zero element \f$x\f$ of the lattice there exists an atom \f$a\f$ such that \f$a\sim x\f$.
    bool atomic = false;
    /// This property is the element 0 of the lattice which satisfies \f$0 \sim x\f$ for all \f$x \ne 0\f$ in the lattice.
    int null = -1;
    /// This property is the element 1 of the lattice which satisfies \f$x \sim 1\f$ for all \f$x \ne 1\f$ in the lattice.
    int unity = -1;

    /// This method calculates which lattice elements are atoms and, once this is done, whether or not the lattice is atomic.
    void compute_atoms();
    /// This method calculates which of the lattice elements are the Lattice::null and Lattice::unity elements.
    void compute_bounds();
    /// This method randomly builds an ordering on the lattice elements such that it satisfies the axioms of a lattice and then calls the Lattice::compute_atoms and Lattice::compute_bounds methods to calculate all the properties of the lattice.
    void initialize();
   public:
    /// The default constructor, which does nothing.
    Lattice();
    /// This constructor accepts as its unique argument the number of elements of the lattice and calls the Poset constructor with this argument, then the Lattice::initialize method. 
    Lattice(int);
    /// The copy constructor that copies over all the class properties.
    Lattice(const Lattice&);
    /// The destructor which does nothing.
    ~Lattice() override;
    /// The overloaded assignment operator that sets all the class properties to those of the source.
    Lattice& operator =(const Lattice&);
    /// This method determines, given the arguments \f$x\f$ and \f$y\f$, the element \f$w\f$ of the lattice such that \f$w \sim x\f$ and \f$w \sim y\f$, while any other element \f$z\f$ in the lattice satisfying these relations also satisfies \f$z \sim w\f$. If the method cannot find such an element it returns N.
    int meet(int,int) const;
    /// This method determines, given the arguments \f$x\f$ and \f$y\f$, the element \f$w\f$ of the lattice such that \f$x \sim w\f$ and \f$y \sim w\f$, while any other element \f$z\f$ in the lattice satisfying these relations also satisfies \f$w \sim z\f$. If the method cannot find such an element it returns N.
    int join(int,int) const;
    /// This method writes the instance properties to a binary disk file and returns the number of bytes written to the file.
    int serialize(std::ofstream&) const override;
    /// This method calls the clear() method on the instance and then reads the properties from a binary disk file and returns the number of bytes read.
    int deserialize(std::ifstream&) override;
    /// This method sets all the class properties, inherited and otherwise, to their default value.
    void clear() override;
    /// This method checks that this object satisfies the fundamental requirements of a lattice, i.e. that every pair of elements the return value of Lattice::meet and Lattice::join is not N.
    bool consistent() const override;
    /// This method calculates the fraction of the potential lattice elements (i.e. thos which are neither null nor atoms themselves) which fail to satisfy the atomic property. 
    double atomicity() const;
  };
}
#endif

