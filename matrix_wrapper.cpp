#include "matrix.cpp"

template class Matrix<int>;
namespace SYNARMOSMA {template unsigned int normalize<int>(Matrix<int>&);}

template class Matrix<NTL::ZZ>;
namespace SYNARMOSMA {template unsigned int normalize<NTL::ZZ>(Matrix<NTL::ZZ>&);}

template class Matrix<double>;

