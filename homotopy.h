#include "group.h"
#include "nexus.h"

#ifndef _homotopyh
#define _homotopyh

namespace SYNARMOSMA {
  /// A class representing the sequence of homotopy groups for a topological space, up to some finite threshold.
  class Homotopy {
   private:
    /// The upper threshold for the homotopy sequence; this value should typically be at 
    /// least two, since for a connected topological space the 0-th homotopy group is the 
    /// trivial group with one element and the 1st homotopy group (the only one which may 
    /// be non-Abelian) is the fundamental group.
    unsigned int ulimit = 0;
    /// The "fitness" of this homotopy sequence, based on an arbitrary formula that measures 
    /// the degree to which the homotopy groups satisfy some structural criteria (e.g. number 
    /// of generators and relations etc.). 
    double fitness = 0.0;
    /// The main property of this class, the actual sequence of Group instances which store the 
    /// homotopy groups; naturally, the length of Homotopy::sequence should be equal to ulimit.
    std::vector<Group> sequence;

    /// This method computes the value of the fitness property, given the formula current implemented: \f$\sum_i \exp[-(n_i-3)^2] + \sin^2(m_i)\f$ where \f$n_i\f$ is the number of generators and \f$m_i\f$ the number of relations for the \f$i\f$-th homotopy group.
    void compute_fitness();
    /// This method selects a random group from the homotopy sequence and alters it, with the severity of the alteration controlled by the value of the method's unique argument which should lie between 0 and 1.
    void mutate(double);
   public:
    /// The default constructor which does nothing.
    Homotopy();
    /// This constructor accepts the argument as the value of Homotopy::ulimit and assembles a set of random instances of the Group class, though the first group is set to be the trivial group and the constructor verifies that the homotopy groups for \f$n>1\f$ are Abelian.
    Homotopy(unsigned int);
    /// The standard copy constructor that copies over the properties from the source instance.
    Homotopy(const Homotopy&);
    /// The standard overloaded assignment operator that copies over the properties from the source instance.
    Homotopy& operator =(const Homotopy&);
    /// The destructor which does nothing.  
    ~Homotopy();
    /// This method clears the Homotopy::sequence vector and sets the other two properties to their default value.
    void clear();
    /// This method writes the groups of the Homotopy::sequence vector using the Group::compact_form method into a string which is returned. 
    std::string write() const;
    /// This method accepts a simplicial complex as its argument and then computes the first two homotopy groups, assuming the complex is connected (otherwise it does nothing).
    void compute(const Nexus*);
    /// This method writes the edge's properties to a binary disk file and returns the number of bytes written to the file.
    int serialize(std::ofstream&) const;
    /// This method calls the clear() method on the instance and then reads the properties from a binary disk file and returns the number of bytes read.
    int deserialize(std::ifstream&);
    /// This method returns the current value of the private property Homotopy::fitness for this instance of the class.
    double get_fitness() const;
    /// This overloaded operator "mates" two homotopy sequences of the same length by computing a random division point and using the first argument's groups initially and the second argument's groups afterwards to assemble the homotopy sequence for the output instance.
    friend Homotopy operator +(const Homotopy&,const Homotopy&);
    /// This overloaded ostream operator writes out the vector Homotopy::sequence using the Group::compact_form method.
    friend std::ostream& operator <<(std::ostream&,const Homotopy&);
  };

  inline double Homotopy::get_fitness() const 
  {
    return fitness;
  }
}
#endif
