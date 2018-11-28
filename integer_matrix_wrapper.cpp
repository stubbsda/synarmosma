#include "integer_matrix.cpp"

namespace SYNARMOSMA {
  template<>
  const int Integer_Matrix<int>::zero = 0;
  template<>
  const int Integer_Matrix<int>::neg1 = -1;
  template<>
  const int Integer_Matrix<int>::unity = 1; 
  
  template<>
  const NTL::ZZ Integer_Matrix<NTL::ZZ>::zero = NTL::to_ZZ(0);
  template<>
  const NTL::ZZ Integer_Matrix<NTL::ZZ>::neg1 = NTL::to_ZZ(-1);
  template<>
  const NTL::ZZ Integer_Matrix<NTL::ZZ>::unity = NTL::to_ZZ(1); 
}

template class SYNARMOSMA::Integer_Matrix<int>;
namespace SYNARMOSMA {template unsigned int normalize<int>(Integer_Matrix<int>&);}

template class SYNARMOSMA::Integer_Matrix<NTL::ZZ>;
namespace SYNARMOSMA {template unsigned int normalize<NTL::ZZ>(Integer_Matrix<NTL::ZZ>&);}

