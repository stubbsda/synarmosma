#include "integer_matrix.cpp"

namespace SYNARMOSMA {
  template<>
  /// This instantiation is necessary for the 
  /// representation of zero as an ordinary signed 
  /// integer.
  const int Integer_Matrix<int>::zero = 0;
  template<>
  /// This instantiation is necessary for the 
  /// representation of -1 as an ordinary signed 
  /// integer.
  const int Integer_Matrix<int>::neg1 = -1;
  template<>
  /// This instantiation is necessary for the 
  /// representation of 1 as an ordinary signed 
  /// integer.
  const int Integer_Matrix<int>::unity = 1; 
  
  template<>
  /// This instantiation is necessary for the 
  /// representation of zero as a multiprecision 
  /// integer using the NTL.
  const NTL::ZZ Integer_Matrix<NTL::ZZ>::zero = NTL::to_ZZ(0);
  template<>
  /// This instantiation is necessary for the 
  /// representation of -1 as a multiprecision 
  /// integer using the NTL.
  const NTL::ZZ Integer_Matrix<NTL::ZZ>::neg1 = NTL::to_ZZ(-1);
  template<>
  /// This instantiation is necessary for the 
  /// representation of 1 as a multiprecision 
  /// integer using the NTL.
  const NTL::ZZ Integer_Matrix<NTL::ZZ>::unity = NTL::to_ZZ(1); 
}

template class SYNARMOSMA::Integer_Matrix<int>;
namespace SYNARMOSMA {template unsigned int normalize<int>(Integer_Matrix<int>&);}

template class SYNARMOSMA::Integer_Matrix<NTL::ZZ>;
namespace SYNARMOSMA {template unsigned int normalize<NTL::ZZ>(Integer_Matrix<NTL::ZZ>&);}

