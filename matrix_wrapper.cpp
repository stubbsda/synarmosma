#include "matrix.cpp"

namespace SYNARMOSMA {
  template<>
  const int Matrix<int>::zero = 0;
  template<>
  const int Matrix<int>::neg1 = -1;
  template<>
  const int Matrix<int>::unity = 1; 

  template<>
  const NTL::ZZ Matrix<NTL::ZZ>::zero = NTL::to_ZZ(0);
  template<>
  const NTL::ZZ Matrix<NTL::ZZ>::neg1 = NTL::to_ZZ(-1);
  template<>
  const NTL::ZZ Matrix<NTL::ZZ>::unity = NTL::to_ZZ(1); 

  template<>
  const double Matrix<double>::zero = 0.0;
  template<>
  const double Matrix<double>::neg1 = -1.0;
  template<>
  const double Matrix<double>::unity = 1.0; 
}

template class SYNARMOSMA::Matrix<int>;
namespace SYNARMOSMA {template unsigned int normalize<int>(Matrix<int>&);}

template class SYNARMOSMA::Matrix<NTL::ZZ>;
namespace SYNARMOSMA {template unsigned int normalize<NTL::ZZ>(Matrix<NTL::ZZ>&);}

template class SYNARMOSMA::Matrix<double>;
