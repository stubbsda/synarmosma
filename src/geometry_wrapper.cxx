#include "geometry.cxx"

namespace SYNARMOSMA {
  template<class kind>
  const double Geometry<kind>::space_quantum = 1e-8;
}

template class Geometry<INT64>;
template class Geometry<double>;