#include "integer_polynomial.cpp"

namespace SYNARMOSMA {
  template<>
  /// This instantiation is necessary for the 
  /// representation of zero as an ordinary signed 
  /// integer.
  const int Integer_Polynomial<int>::zero = 0;
  template<>
  /// This instantiation is necessary for the 
  /// representation of -1 as an ordinary signed 
  /// integer.
  const int Integer_Polynomial<int>::neg1 = -1;
  template<>
  /// This instantiation is necessary for the 
  /// representation of 1 as an ordinary signed 
  /// integer.
  const int Integer_Polynomial<int>::unity = 1; 
  
  template<>
  /// This instantiation is necessary for the 
  /// representation of zero as a multiprecision 
  /// integer using the NTL.
  const NTL::ZZ Integer_Polynomial<NTL::ZZ>::zero = NTL::to_ZZ(0);
  template<>
  /// This instantiation is necessary for the 
  /// representation of -1 as a multiprecision 
  /// integer using the NTL.
  const NTL::ZZ Integer_Polynomial<NTL::ZZ>::neg1 = NTL::to_ZZ(-1);
  template<>
  /// This instantiation is necessary for the 
  /// representation of 1 as a multiprecision 
  /// integer using the NTL.
  const NTL::ZZ Integer_Polynomial<NTL::ZZ>::unity = NTL::to_ZZ(1); 
}

template class SYNARMOSMA::Integer_Polynomial<int>;
template class SYNARMOSMA::Integer_Polynomial<NTL::ZZ>;

