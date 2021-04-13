#include "global.h"

#ifndef _vertexh
#define _vertexh

namespace SYNARMOSMA {
  template<class kind>
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
    kind energy = kind(0);
    /// This integer set contains the index of the vertices connected to this vertex by 
    /// outgoing directed edges. To the extent that the direction of an edge has a causal 
    /// sense, this set thus represents the immediate future of this vertex-event. 
    std::set<int> posterior;
    /// This integer set contains the index of the vertices connected to this vertex by 
    /// incoming directed edges. To the extent that the direction of an edge has a causal 
    /// sense, this set thus represents the immediate past of this vertex-event. 
    std::set<int> anterior;
    /// This integer set contains the index of all vertices connected to this vertex, thus the 
    /// elements of the Vertex::anterior and Vertex::posterior properties but also vertices that 
    /// are connected by an undirected (Relation::disparate) edge.
    std::set<int> neighbours;
    /// This property is a set of the index of all of the edges (as stored by the
    /// super class) which contain this vertex as an endpoint.
    std::set<int> entourage;
    /// This constant represents the smallest possible energy value for a 
    /// vertex; it is only meaningful when this template class is instantiated 
    /// with a discrete base type.
    static const double energy_quantum;

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
    /// This method sets the energy value using a double, which if the base type is UINT64 is then divided by the energy quantum to get an unsigned 64 bit integer.
    void set_energy(double);
    /// This method increments the energy value using a double, which if the base type is UINT64 is then divided by the energy quantum to get an unsigned 64 bit integer that can be added to the current value.
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

  template<>
  /// This method has been specialized for the UINT64 base type, since in this case the method can verify exactly that the energy is zero.
  inline bool Vertex<UINT64>::zero_energy() const
  {
    return (energy == 0);
  }

  template<class kind>
  inline bool Vertex<kind>::zero_energy() const
  {
    bool output = (std::abs(energy) < std::numeric_limits<double>::epsilon()) ? true : false;
    return output;
  }

  template<class kind>
  inline void Vertex<kind>::nullify_energy()
  {
    energy = kind(0);
  }

  template<>
  /// This method has been specialized for the UINT64 base type, since in this case the argument must be divided by Vertex::energy_quantum and the fractional part truncated before incrementing Vertex::energy by it.
  inline void Vertex<UINT64>::increment_energy(double E)
  {
    if (E <= std::numeric_limits<double>::epsilon()) throw std::invalid_argument("Vertex::increment_energy argument must be greater than zero!");
    energy += UINT64(E/energy_quantum);
  }

  template<class kind>
  inline void Vertex<kind>::increment_energy(double E)
  {
    if (E <= std::numeric_limits<double>::epsilon()) throw std::invalid_argument("Vertex::increment_energy argument must be greater than zero!");
    energy += E;
  }

  template<>
  /// This method has been specialized for the UINT64 base type, since in this case the argument must be divided by Vertex::energy_quantum and the fractional part truncated before setting Vertex::energy to it.
  inline void Vertex<UINT64>::set_energy(double E)
  {
    if (E <= std::numeric_limits<double>::epsilon()) throw std::invalid_argument("Vertex::set_energy argument must be greater than zero!");
    energy = UINT64(E/energy_quantum);
  }

  template<class kind>
  inline void Vertex<kind>::set_energy(double E)
  {
    if (E <= std::numeric_limits<double>::epsilon()) throw std::invalid_argument("Vertex::set_energy argument must be greater than zero!");
    energy = E;
  }

  template<>
  /// This method has been specialized for the UINT64 base type, since in this case the Vertex::energy value must be converted to a double and then multiplied by Vertex::energy_quantum prior to be being returned.
  inline double Vertex<UINT64>::get_energy() const
  {
    return energy_quantum*double(energy);
  }

  template<class kind>
  inline double Vertex<kind>::get_energy() const
  {
    return energy;
  }
}
#endif 

