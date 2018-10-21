#include "global.h"

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
    /// number, represents the energy of a vertex and should therefore always be positive. 
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
    Vertex();
    Vertex(const Vertex&);
    virtual ~Vertex();
    Vertex& operator =(const Vertex&);
    inline bool zero_energy() const;
    inline double get_energy() const;
    inline void set_energy(double);
    inline void increment_energy(double);
    inline void nullify_energy();
    virtual int serialize(std::ofstream&) const;
    virtual int deserialize(std::ifstream&);
    virtual void clear();
  };

  bool Vertex::zero_energy() const
  {
#ifdef DISCRETE
    bool output = (energy == 0) ? true : false;
#else
    bool output = (std::abs(energy) < std::numeric_limits<double>::epsilon()) ? true : false;
#endif
    return output;
  }

  void Vertex::nullify_energy()
  {
#ifdef DISCRETE
    energy = 0;
#else
    energy = 0.0;
#endif
  }

  void Vertex::increment_energy(double E)
  {
#ifdef DEBUG
    assert(E > std::numeric_limits<double>::epsilon());
#endif
#ifdef DISCRETE
    energy += UINT64(E/energy_quantum);
#else
    energy += E;
#endif
  }

  void Vertex::set_energy(double E)
  {
#ifdef DEBUG
    assert(E > std::numeric_limits<double>::epsilon());
#endif
#ifdef DISCRETE
    energy = UINT64(E/energy_quantum);
#else
    energy = E;
#endif
  }

  double Vertex::get_energy() const
  {
#ifdef DISCRETE
    return energy_quantum*double(energy);
#else
    return energy;
#endif
  }
}
#endif 

