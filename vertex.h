#include "random.h"

#ifndef _vertexh
#define _vertexh

namespace SYNARMOSMA {
  /// A class representing the vertex, node or 0-simplex of a simplicial complex of dimension greater than zero or the elements of a set.
  class Vertex {
   protected:
    /// This integer property encodes the relaxation or time step at which the vertex was 
    /// created. 
    int incept = -1;
    /// This integer property reflects the highest dimensional simplex to which this vertex 
    /// belongs. 
    int topological_dimension = 0;
    /// This property, either an unsigned 64 bit integer or a double precision floating point 
    /// number, represents the energy of a vertex and should therefore always be non-negative. 
#ifdef DISCRETE
    UINT64 energy = 0;
#else
    double energy = 0.0;
#endif
    /// This integer set contains the index of the vertices connected to this vertex by 
    /// outgoing directed edges. To the extent that the direction of an edge has a causal 
    /// sense, this set thus represents the immediate future of this vertex-event. 
    std::set<int> posterior;
    /// This integer set contains the index of the vertices connected to this vertex by 
    /// incoming directed edges. To the extent that the direction of an edge has a causal 
    /// sense, this set thus represents the immediate past of this vertex-event. 
    std::set<int> anterior;
    /// This integer set contains the index of all vertices connected to this vertex, thus the 
    /// elements of the "anterior" and "posterior" properties but also vertices connected by an 
    /// undirected (spacelike) edge.
    std::set<int> neighbours;
    /// This property is a set of the index of all of the edges (as stored by the
    /// super class) which contain this vertex as an endpoint.
    std::set<int> entourage;

   public:
    /// The default constructor that sets all properties to their default value.
    Vertex();
    /// Standard copy constructor for this class.
    Vertex(const Vertex&);
    /// Standard destructor for this class, virtual since the class is expected to have sub-classes.
    virtual ~Vertex();
    /// The assignment operator that copies over all class properties from the source to the target.
    Vertex& operator =(const Vertex&);
    /// This method tests if the energy property of this instance is zero.
    bool zero_energy() const;
    /// This returns the energy of the vertex as a double precision number regardless of the storage format (UINT64 or double).
    double get_energy() const;
    /// This method sets the energy value using a double, which if DISCRETE is then divided by the energy quantum to get an unsigned 64 bit integer.
    void set_energy(double);
    /// This method increments the energy value using a double, which if DISCRETE is then divided by the energy quantum to get an unsigned 64 bit integer that can be added to the current value.
    void increment_energy(double);
    /// This method sets the energy property equal to zero.
    void nullify_energy();
    /// This method writes the instance properties to a binary disk file and returns the number of bytes written to the file.
    virtual int serialize(std::ofstream&) const;
    /// This method calls the clear() method on the instance and then reads the properties from a binary disk file and returns the number of bytes read.
    virtual int deserialize(std::ifstream&);
    /// This method restores all the properties of the instance to their initialization value.
    virtual void clear();
  };

  inline bool Vertex::zero_energy() const
  {
#ifdef DISCRETE
    bool output = (energy == 0) ? true : false;
#else
    bool output = (std::abs(energy) < std::numeric_limits<double>::epsilon()) ? true : false;
#endif
    return output;
  }

  inline void Vertex::nullify_energy()
  {
#ifdef DISCRETE
    energy = 0;
#else
    energy = 0.0;
#endif
  }

  inline void Vertex::increment_energy(double E)
  {
    if (E <= std::numeric_limits<double>::epsilon()) throw std::invalid_argument("Vertex::increment_energy argument must be greater than zero!");
#ifdef DISCRETE
    energy += UINT64(E/energy_quantum);
#else
    energy += E;
#endif
  }

  inline void Vertex::set_energy(double E)
  {
    if (E <= std::numeric_limits<double>::epsilon()) throw std::invalid_argument("Vertex::set_energy argument must be greater than zero!");
#ifdef DISCRETE
    energy = UINT64(E/energy_quantum);
#else
    energy = E;
#endif
  }

  inline double Vertex::get_energy() const
  {
#ifdef DISCRETE
    return energy_quantum*double(energy);
#else
    return energy;
#endif
  }
}
#endif 

