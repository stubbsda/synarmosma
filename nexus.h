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
    Nexus(int);
    Nexus(const Nexus&);
    Nexus& operator =(const Nexus&);
    /// The destructor which frees the memory of the Nexus::elements and Nexus::index_table properties if Nexus::dimension is greater than -1.
    ~Nexus() override;
    int serialize(std::ofstream&) const override;
    int deserialize(std::ifstream&) override;
    bool orientable() const;
    bool pseudomanifold(bool*) const;
    void surface_construction(std::string&);
    void clear() override;
    void initialize(int);
    void initialize(int,int);
    inline bool paste(const std::set<int>& vx) {return paste(Cell(vx));};
    /// This method adds a new Cell to this Nexus instance and returns true if it is genuinely new, false otherwise. Note that this method does not call regularization(), so it is the caller's responsibility to ensure that the Nexus is explicitly regularized.
    bool paste(const Cell&);
    /// This method ensures that this Nexus instance satisfies the entailment property for an abstract simplicial complex, i.e. if a D-simplex is part of the complex then so are all of its sub-simplices. The method returns the number of simplices that had to be added to regularize the complex.
    int regularization();
    void ascend(int,int,std::vector<Cell>&) const;
    void star(const std::set<std::set<int> >&,std::vector<Cell>*) const;
    void link(const std::set<std::set<int> >&,std::vector<Cell>*) const;
    void closure(const std::set<std::set<int> >&,Nexus*,int*) const;
    void compute_entourages();
    void compute_neighbours();
    /// This method returns the dimensionality of this instance of the Nexus class.  
    inline int get_dimension() const {return dimension;};
    /// This method returns the number of simplices of a given dimension in this instance. 
    inline unsigned int get_length(int) const;
    /// This method returns the index of a Cell of a given dimension in the relevant Nexus::elements vector and -1 if this Cell can't be found.
    inline int get_index(std::set<int>&) const;
    inline void get_elements(int D,int n,int* vx) const {elements[D][n].get_vertices(vx);};
    inline void get_elements(int D,int n,std::set<int>& vx) const {elements[D][n].get_vertices(vx);};
    inline void get_entourage(int D,int n,std::set<int>& vx) const {elements[D][n].get_entourage(vx);};
  };

  unsigned int Nexus::get_length(int D) const 
  {
    if (D < 1) throw std::invalid_argument("The dimensionality in Nexus::get_length must be greater than zero!"); 
    return elements[D].size();
  }

  int Nexus::get_index(std::set<int>& S) const 
  {
    unsigned int D = S.size() - 1; 
    if (D < 1) throw std::invalid_argument("The dimensionality of the simplex in Nexus::get_index must be greater than zero!"); 
    hash_map::const_iterator qt = index_table[D].find(S); 
    if (qt != index_table[D].end()) return qt->second;
    return -1;
  }
}
#endif


