#include "multitime.cxx"

namespace SYNARMOSMA {
  template<class kind>
  const int Multitime<kind>::tdimension;

  template<class kind>
  const double Multitime<kind>::time_quantum = 1e-8;
}

template class Multitime<INT64>;
template class Multitime<double>;