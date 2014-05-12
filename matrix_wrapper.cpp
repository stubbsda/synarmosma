#include "matrix.cpp"

template class Matrix<int>;
template unsigned int normalize<int>(Matrix<int>&);

template class Matrix<NTL::ZZ>;
template unsigned int normalize<NTL::ZZ>(Matrix<NTL::ZZ>&);

template class Matrix<double>;

