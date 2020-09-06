#include "variety.cxx"

namespace SYNARMOSMA {
  template<>
  /// This instantiation is necessary for the 
  /// representation of zero as a rational number, 
  /// i.e. 0/1.
  const Rational Variety<Rational>::zero = Rational(0);

  template<>
  /// This instantiation is necessary for the 
  /// representation of zero as an ordinary unsigned  
  /// integer.
  const unsigned int Variety<unsigned int>::zero = 0;

  template<>
  /// This instantiation is necessary for the 
  /// representation of zero as an ordinary signed  
  /// integer.
  const int Variety<int>::zero = 0;

  template<>
  /// This instantiation is necessary for the 
  /// representation of zero as a multiprecision 
  /// integer using the NTL.
  const NTL::ZZ Variety<NTL::ZZ>::zero = NTL::to_ZZ(0);
}

template class SYNARMOSMA::Variety<Rational>;
template class SYNARMOSMA::Variety<unsigned int>;
template class SYNARMOSMA::Variety<int>;
template class SYNARMOSMA::Variety<NTL::ZZ>;

