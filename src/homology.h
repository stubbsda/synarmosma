#include "nexus.h"
#include "integer_matrix.h"

#ifndef _homologyh
#define _homologyh

namespace SYNARMOSMA {
  /// A class representing the finite sequence of homology groups for a topological space, a sequence whose length is equal to the space's dimension. 
  class Homology {
   public:
    /// This enumerated class lists the three kind of domains over which the homology can be computed: native 32 bit integers, multiprecision 
    /// integers (NTL::ZZ) and mod2. 
    enum class Field
    {
        int32,
        multiprecision,
        mod2
    };
    /// This enumerated class lists the two methods available for computing the homology, either using the computer algebra software GAP (if it 
    /// is available) or a native technique in the Homology class itself.
    enum class Method
    {
        gap,
        native
    };

   private:
    /// The Betti numbers for the topological space - the homology groups are all finitely-generated 
    /// Abelian groups and so we can store them by noting the rank of each group (i.e. the Betti number) 
    /// and the vector of torsions. 
    std::vector<unsigned int> betti_number;
    /// The torsion vector for each homology group, as the torsion for a finitely-generated Abelian 
    /// group can be written as \f$T = \mathbf{Z}_{n_1} \oplus \mathbf{Z}_{n_2} \oplus \dots \oplus 
    /// \mathbf{Z}_{n_k}\f$ so we need only store the values of \f$(n_1,n_2,\dots,n_k)\f$ to know the 
    /// torsion, if any.
    std::vector<std::vector<unsigned int> > torsion;
    /// This property determines which of the two potential methods should be used to compute the homology 
    /// groups - the GAP software or by the built-in solver for calculating the Jordan normal form of the 
    /// integer matrices.
    Method method = Method::native;
    /// This property controls the domain over which the homology is calculated - either the integers, 
    /// represented using the native "int" type or the multiprecision NTL::ZZ type, or else the Galois 
    /// field GF(2).
    Field field = Field::int32;

    /// This method computes the integral (so Homology::field equals int32 or multiprecision) homology of an instance of the Nexus class when Homology::method equals native.
    void compute_integral_native(const Nexus*);
    /// This method computes the integral or mod2 homology of an instance of the Nexus class without using GAP; the simpler case of Homology::field equals mod2 is handled by this method itself while for the integral homology the compute_integral_native() method is called. 
    void compute_native(const Nexus*);
    /// This method computes the integral or mod2 homology of an instance of the Nexus class using the GAP software; in particular, for integral homology the GAP package "simpcomp" is required.
    void compute_gap(const Nexus*);
   public:
    /// The default constructor, which does nothing.
    Homology();
    /// The standard constructor whose two arguments set the values for the Homology::field and Homology::method properties.
    Homology(Field,Method);
    /// The copy constructor that copies over all the class properties.
    Homology(const Homology&);
    /// The overloaded assignment operator that sets all the class properties to those of the source.
    Homology& operator =(const Homology&);
    /// The destructor which does nothing.
    ~Homology();
    /// This method writes out the sequence of homology groups as a string, using a comma-separated list where each group has the format \f$[r; t_1,t_2,\dots,t_k]\f$ where \f$r\f$ is the Betti number and the \f$t_i\f$ are the elements of the torsion vector.
    std::string write() const;
    /// The standard "set" method that gives a value to the Homology::method property.
    void set_method(Method);
    /// The standard "set" method that gives a value to the Homology::field property.
    void set_field(Field);
    /// The standard "get" method that returns the current value of the Homology::method property.
    Method get_method() const;
    /// The standard "get" method that returns the current value of the Homology::field property.
    Field get_field() const;
    /// This method sets the argument to the vector of Betti numbers that have been calculated.
    void get_betti_numbers(std::vector<unsigned int>&) const; 
    /// This method clears the Homology::betti_number and Homology::torsion vectors and sets the other class properties to their default value.
    void clear();
    /// This method computes the homology groups of an instance of the Nexus class, i.e. an abstract simplicial complex.
    void compute(const Nexus*);
    /// This method writes the instance properties to a binary disk file and returns the number of bytes written to the file.
    int serialize(std::ofstream&) const;
    /// This method calls the clear() method on the instance and then reads the properties from a binary disk file and returns the number of bytes read.
    int deserialize(std::ifstream&);
  };

  inline void Homology::set_method(Method m) 
  {
    method = m;
  }

  inline Homology::Method Homology::get_method() const 
  {
    return method;
  }

  inline void Homology::set_field(Field f) 
  {
    field = f;
  }

  inline Homology::Field Homology::get_field() const 
  {
    return field;
  }

  inline void Homology::get_betti_numbers(std::vector<unsigned int>& output) const
  {
    output = betti_number;
  }
}
#endif 
