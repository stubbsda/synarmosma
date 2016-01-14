/*
  Copyright 2014 Daniel Stubbs

  This file is part of Synarmosma.

  Synarmosma is free software: you can redistribute it and/or modify 
  it under the terms of the GNU General Public License as published by 
  the Free Software Foundation, either version 3 of the License, or 
  (at your option) any later version.

  Synarmosma is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with Synarmosma.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "proposition.h"

#ifndef _vertexh
#define _vertexh

namespace SYNARMOSMA {
  class Vertex {
   protected:
    int incept;
    int topological_dimension;
#ifdef DISCRETE
    UINT64 energy;
#else
    double energy;
#endif
    Proposition theorem;
    std::set<int> posterior,anterior;  
    std::set<int> neighbours;
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
    virtual void serialize(std::ofstream&) const;
    virtual void deserialize(std::ifstream&);
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

