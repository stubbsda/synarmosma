#include "variety.cpp"

namespace SYNARMOSMA {
  template<>
  const Rational Variety<Rational>::zero = Rational(0);

  template<>
  const unsigned int Variety<unsigned int>::zero = 0;

  template<>
  const int Variety<int>::zero = 0;

  template<>
  const NTL::ZZ Variety<NTL::ZZ>::zero = NTL::to_ZZ(0);
}

template class SYNARMOSMA::Variety<Rational>;
template class SYNARMOSMA::Variety<unsigned int>;
template class SYNARMOSMA::Variety<int>;
template class SYNARMOSMA::Variety<NTL::ZZ>;

