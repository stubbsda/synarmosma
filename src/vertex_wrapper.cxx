#include "vertex.cxx"

namespace SYNARMOSMA {
  template<class kind>
  const double Vertex<kind>::energy_quantum = 1e-7;
}

template class SYNARMOSMA::Vertex<UINT64>;
template class SYNARMOSMA::Vertex<double>;
