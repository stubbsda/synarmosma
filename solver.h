#include <complex>
#include <vector>
#include "synarmosma/matrix.h"

class Solver {
 protected:
  unsigned int max_its;
  unsigned int dim;
  unsigned int nnonzero;
  double epsilon;  
  double t;
  std::vector<std::complex<double> > c_solution;
  SYNARMOSMA::Matrix<std::complex<double> >* J;

  bool forward_step();
  void set_default_values();
  void compute_jacobian(const std::vector<std::complex<double> >&);
  void linear_solver(const std::vector<std::complex<double> >&,const std::vector<std::complex<double> >&,std::vector<std::complex<double> >&) const;
  std::complex<double> a1(double) const;
  std::complex<double> a2(double) const;
  virtual void F(const std::vector<std::complex<double> >&,std::vector<std::complex<double> >&) const = 0;
 public:
  Solver(int);
  ~Solver();
  bool solve(std::vector<std::complex<double> >&);
};

