#include "polynomial.cpp"

template class Polynomial<Rational>;
template class Polynomial<unsigned int>;
template class Polynomial<signed int>;
template class Polynomial<double>;
template class Polynomial<std::complex<double> >;
template class Polynomial<NTL::ZZ>;

