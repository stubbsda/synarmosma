#include "cell.h"
#include "graph.h"

#ifndef _nexush
#define _nexush

namespace SYNARMOSMA {
  /// A class representing an abstract simplicial complex, using the Cell class to store higher-dimensional elements.
  class Nexus: public Schema {
   protected:
    /// This integer property stores the dimension of the abstract simplicial complex, i.e. 
    /// the length of the array Nexus::elements less one; it is set when the instance of the 
    /// Nexus is constructed and cannot be changed. 
    int dimension = -1;
    /// This property is an array of vectors of Cell type, each such vector corresponding to a 
    /// dimension. The 0-dimensional vector, Nexus::elements[0], is always empty as the number 
    /// of vertices is known from the Schema::nvertex property that this class inherits.  
    std::vector<Cell>* elements;
    /// This property stores the mapping between the set of vertices of a d-dimensional Cell and its 
    /// index in the vector Nexus::elements[d].
    hash_map* index_table;

   public:
    /// The default constructor which does nothing.
    Nexus();
    /// The standard constructor for the Nexus class, in which the argument stores the value of the dimension which must be greater than zero.
    Nexus(int);
    /// The standard copy constructor - it calls the clear() method and then copies over all properties from the source instance to this one.
    Nexus(const Nexus&);
    /// The overloaded assignment operator for this class, which first calls clear() and then behaves exactly like the copy constructor for this class.
    Nexus& operator =(const Nexus&);
    /// The destructor which frees the memory of the Nexus::elements and Nexus::index_table properties if Nexus::dimension is greater than -1.
    ~Nexus() override;
    /// This method writes the instance properties to a binary disk file and returns the number of bytes written to the file.
    int serialize(std::ofstream&) const override;
    /// This method calls the clear() method on the instance and then reads the properties from a binary disk file and returns the number of bytes read.
    int deserialize(std::ifstream&) override;
    /// This method determines if the abstract simplicial complex, assumed to be pseudomanifold, is orientable and returns true if this is so.
    bool orientable() const;
    /// This method determines if the abstract simplicial complex represented by this Nexus instance is a pseudomanifold and returns true if this is so; the method's argument is a Boolean which is set to true if the Nexus instance is a pseudomanifold-with-boundary.
    bool pseudomanifold(bool*) const;
    /// This method builds a minimal simplicial complex reflecting the topology of one of a fixed set of know surface topologies: SPHERE (4), PROJECTIVE_PLANE (6), TORUS (9) and MÖBIUS_STRIP (8); the number in brackets indicates the number of vertices needed to represent this surface. 
    void surface_construction(std::string&);
    /// This method restores the inherited Schema properties to their default value and, if Nexus::dimension is greater than -1, frees the memory for Nexus::elements and Nexus::index_table, then sets Nexus::dimension to -1. 
    void clear() override;
    /// This method calls clear(), sets the value of Nexus::dimension to the argument and allocates the Nexus::elements and Nexus::index_table properties. 
    void initialize(int);
    /// This method calls clear(), sets the value of Nexus::dimension to the first argument and the inherited nvertex property to the second argument. The method then allocates the Schema and Nexus properties.
    void initialize(int,int);
    /// This method functions like the other paste() method but uses a set containing the vertices that should be added to the Nexus instance.
    inline bool paste(const std::set<int>& vx) {return paste(Cell(vx));};
    /// This method adds a new Cell to this Nexus instance and returns true if it is genuinely new, false otherwise. Note that this method does not call regularization(), so it is the caller's responsibility to ensure that the Nexus is explicitly regularized.
    bool paste(const Cell&);
    /// This method ensures that this Nexus instance satisfies the entailment property for an abstract simplicial complex, i.e. if a D-simplex is part of the complex then so are all of its sub-simplices. The method returns the number of simplices that had to be added to regularize the complex.
    int regularization();
    /// This method computes via recursion the total entourage of a given d-simplex S, i.e. all n-simplices with d < n <= dimension that contain S. The first argument is the dimension d, the second the index of the simplex in Nexus::elements[d] and the final argument is the vector that will stores all of these instances of the Cell class. 
    void ascend(int,int,std::vector<Cell>&) const;
    /// This method computes the star of a set of simplices, contained in the first argument as a set of vertex sets. The star of the set S is the union of the stars of each simplex in S. For a single simplex s, the star of s is the set of simplices having s as a face; note that the star of S is generally not a simplicial complex itself. For this reason the output - the method's second argument - is stored as simply an array of Cell vectors, one for each dimension. 
    void star(const std::set<std::set<int> >&,std::vector<Cell>*) const;
    /// This method computes the link of a set of simplices, contained in the first argument as a set of vertex sets. The link of of a set S equals closure(star(S)) − star(closure(S)); it is the closed star of S minus the stars of all faces of S. Like the star operator, the output of the link operator generally isn't a simplicial complex and so the output, the method's second argument, is stored as an array of Cell vectors, one for each dimension. 
    void link(const std::set<std::set<int> >&,std::vector<Cell>*) const;
    /// This method computes the closure of a set of simplices, contained in the first argument as a set of vertex sets. The closure of the set S is the smallest simplicial subcomplex of this Nexus instance that contains each simplex in S; it is obtained by repeatedly adding to S each face of every simplex in S. The second argument is the closure itself and the final argument is the offset vector for the vertices of this new Nexus instance and this one.
    void closure(const std::set<std::set<int> >&,Nexus*,int*) const;
    /// This method verifies that the simplicial complex is consistent, i.e. that it satisfies the entailment property, that Nexus::elements[0] and Nexus::index_table are empty and that as a Schema it is consistent.
    virtual bool consistent() const override;
    /// This method computes the set of entourages for the Nexus, i.e. verifies that the given sub-simplex exists and records its index value. 
    void compute_entourages();
    /// This method clears the neighbour sets and then recomputes them on the basis of the set of 1-simplices contained in the Nexus instance.
    void compute_neighbours();
    /// This method returns the dimensionality of this instance of the Nexus class.  
    inline int get_dimension() const {return dimension;};
    /// This method returns the number of simplices of a given dimension in this instance. 
    inline unsigned int get_length(int) const;
    /// This method returns the index of a Cell of a given dimension in the relevant Nexus::elements vector and -1 if this Cell can't be found.
    inline int get_index(std::set<int>&) const;
    /// This method puts the vertices belonging to the Cell specified by the first two arguments - the dimension d and index in Nexus::elements[d] - into the integer array in the final argument.
    inline void get_elements(int,int,int*) const;
    /// This method puts the vertices belonging to the Cell specified by the first two arguments - the dimension d and index in Nexus::elements[d] - into the integer set in the final argument.
    inline void get_elements(int,int,std::set<int>&) const;
    /// This method puts the entourage associated with the Cell specified by the first two arguments - the dimension d and index in Nexus::elements[d] - into the integer set in the final argument.
    inline void get_entourage(int,int,std::set<int>&) const;
  };

  unsigned int Nexus::get_length(int D) const 
  {
    if (D < 1 || D > dimension) throw std::invalid_argument("Illegal dimension value in Nexus::get_length!"); 
    return elements[D].size();
  }

  int Nexus::get_index(std::set<int>& S) const 
  {
    unsigned int D = S.size() - 1; 
    if (D < 1 || D > (unsigned) dimension) throw std::invalid_argument("Illegal dimension value in Nexus::get_index!"); 
    hash_map::const_iterator qt = index_table[D].find(S); 
    if (qt != index_table[D].end()) return qt->second;
    return -1;
  }

  void Nexus::get_elements(int D,int n,int* vx) const 
  {
    if (D < 1 || D > dimension) throw std::invalid_argument("Illegal dimension value in Nexus::get_elements!"); 
    if (n < 0 || n >= (signed) elements[D].size()) throw std::invalid_argument("The simplex specified in Nexus::get_elements does not exist!");

    elements[D][n].get_vertices(vx);
  }

  void Nexus::get_elements(int D,int n,std::set<int>& vx) const 
  {
    if (D < 1 || D > dimension) throw std::invalid_argument("Illegal dimension value in Nexus::get_elements!"); 
    if (n < 0 || n >= (signed) elements[D].size()) throw std::invalid_argument("The simplex specified in Nexus::get_elements does not exist!");

    elements[D][n].get_vertices(vx);
  }

  void Nexus::get_entourage(int D,int n,std::set<int>& vx) const 
  {
    if (D < 1 || D > dimension) throw std::invalid_argument("Illegal dimension value in Nexus::get_entourage!"); 
    if (n < 0 || n >= (signed) elements[D].size()) throw std::invalid_argument("The simplex specified in Nexus::get_entourage does not exist!");

    elements[D][n].get_entourage(vx);
  }
}
#endif


