#include "word.h"

#ifndef _grouph
#define _grouph

namespace SYNARMOSMA {
  /// A class representing a combinatorial presentation of a group, using the library's Word class to store the relations. 
  class Group {
   protected:
    /// This property is the number of generators for this presentation of 
    /// the group.
    unsigned int ngenerator = 0;
    /// This property is the number of distinct elements of the group; it is 
    /// only meaningful when the finite property is true.
    unsigned int cardinality = 0;
    /// This property is the rank of the group and only makes sense in the 
    /// context of an Abelian group, for which we have a canonical structure 
    /// based on the rank and torsion. If the rank > 0, the group is not finite.
    unsigned int rank = 0;
    /// This property reflects whether or not all the elements of the group 
    /// commute with one another for the group operation, \f$a\cdot b = b\cdot a\f$; 
    /// when it is false the rank, torsion and solvable properties are meaningless.
    bool abelian = false;
    /// This property reflects whether or not the group has a finite number of 
    /// elements and so whether or not the cardinality property is meaningful.
    bool finite = false;
    /// This property is relevant only for Abelian groups and is true when the 
    /// group can be written as a subnormal series whose quotient groups are all 
    /// Abelian. 
    bool solvable = false;
    /// This property indicates whether or not this group presentation has any 
    /// relations; the group is free when there are none. 
    bool free = true;
    /// This is true when this group belongs to the class of braid groups, which 
    /// means that it has a solvable word problem unlike the general case.
    bool braid = false;
    /// Like the property "rank", the torsion is only meaningful in the context of 
    /// an Abelian group and reflects the extent to which the group possesses a finite 
    /// cyclic structure. If an Abelian group has a rank = 0 and an empty torsion 
    /// vector, then it is the trivial group. 
    std::vector<unsigned int> torsion;
    /// The relations for the group presentation, i.e. a set of words on the group's 
    /// generators, with the understanding that each word is equal to the identity element.
    std::vector<Word> relations;

    /// This method determines the rank of the group, first checking that it is Abelian; it currently sets the rank to zero if the group is finite and to the number of generators otherwise, though this is just an upper bound for the rank.
    void compute_rank();
    /// This method makes a new Abelian group from the current instance by adding the relation \f$aba^{-1}b^{-1}\f$ for each distinct pair of generators \f$a\f$ and \f$b\f$.
    Group abelianize() const;
   public:
    /// The default constructor which simply initializes the properties with their default values.
    Group();
    /// This constructor accepts as its argument the number of generators and builds the free group on n generators; if n = 0 this is the trivial group, if n = 1 the infinite cyclic group of integers under addition. 
    Group(unsigned int);
    /// This constructor builds the group with a given number of generators (first argument) and a set of relations (second argument).
    Group(unsigned int,const std::vector<Word>&);
    /// This constructor accepts two arguments, the number of generators and the number of relations, and constructs a random group presentation obeying these constraints.
    Group(unsigned int,unsigned int);
    /// This is the constructor for an Abelian group in canonical form, with the first argument the rank and the second the torsion coefficients.
    Group(unsigned int,const std::vector<unsigned int>&);
    /// This constructor accepts a series of potential values for the named group type (first argument) as well as an index, permitting the creation of a canonical presentation for such groups as the dihedral or cyclic.
    Group(std::string&,unsigned int);
    /// The standard copy constructor.
    Group(const Group&);
    /// The standard deconstructor that does nothing in this case.
    ~Group();
    /// The overloaded assignment operator that copies over all the group properties.
    Group& operator =(const Group&);
    /// This method simply returns the value of the rank property.
    inline unsigned int get_rank() const {return rank;};
    /// This method computes all of the letters in the collection of group relations (the first argument) and returns true if there are no letters which are not generators, false otherwise.
    bool consistent(std::set<unsigned int>&) const;
    /// This method initializes the group presentation after the value of ngenerator has been specified (in the constructor) and where the argument is the number of relations, which are constructed randomly.
    void initialize(unsigned int);
    /// This method initializes the group presentation, where the first argument is the number of generators and the second the accompanying relations.
    void initialize(unsigned int,const std::vector<Word>&);
    /// This method initializes the group structure for an Abelian group, where the first argument is the rank and the second the torsion coefficients.
    void initialize(unsigned int,const std::vector<unsigned int>&);
    /// This method tests if the two words are equivalent in this group, returning true if so; note that using this method is only meaningful if this is a braid group. 
    bool equivalent(const Word&,const Word&) const;
    /// This method attempts to simplify the group presentation by reducing the words in its relations to their most elementary form, such as eliminating duplicate relations, normalizing the words and using a relation like \f$ab^{-1}\f$ (which implies \f$a=b\f$) to reduce the number of generators.
    void reduce();
    /// This method returns all the properties to their default values and empties the "relations" vector.
    void clear();
    /// This method writes the instance properties to a binary disk file and returns the number of bytes written to the file.
    int serialize(std::ofstream&) const;
    /// This method calls the clear() method on the instance and then reads the properties from a binary disk file and returns the number of bytes read.
    int deserialize(std::ifstream&);
    /// This method writes out the group presentation in the standard mathematical notation \f$\langle \{x_1,\dots,x_n\} \mid \{r_1,\dots,r_m\}\rangle\f$ for non-Abelian groups, while Abelian ones are written out using the rank and torsion coefficients.
    std::string compact_form() const;
    /// The overloaded ostream operator writes out the group structure, beginning with a line listing the scalar properties (Abelian, free, solvable etc.) followed by a second line with the group's combinatorial presentation.
    friend std::ostream& operator <<(std::ostream&,const Group&);
    friend class Homotopy;
  };
}
#endif
