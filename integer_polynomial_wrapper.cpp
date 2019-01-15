#include "integer_polynomial.cpp"

namespace SYNARMOSMA {
  template<>
  const int Integer_Polynomial<int>::zero = 0;
  template<>
  const int Integer_Polynomial<int>::neg1 = -1;
  template<>
  const int Integer_Polynomial<int>::unity = 1; 
  
  template<>
  const NTL::ZZ Integer_Polynomial<NTL::ZZ>::zero = NTL::to_ZZ(0);
  template<>
  const NTL::ZZ Integer_Polynomial<NTL::ZZ>::neg1 = NTL::to_ZZ(-1);
  template<>
  const NTL::ZZ Integer_Polynomial<NTL::ZZ>::unity = NTL::to_ZZ(1); 
}

template class SYNARMOSMA::Integer_Polynomial<int>;
template class SYNARMOSMA::Integer_Polynomial<NTL::ZZ>;

