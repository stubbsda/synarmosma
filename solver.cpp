#include "solver.h"

using namespace SYNARMOSMA;

extern Random RND;

template<class kind>
Solver<kind>::Solver(int N)
{
  set_default_values();
  dim = N;
  J = new SYNARMOSMA::Matrix<kind>(dim);
}

template<class kind>
Solver<kind>::~Solver()
{
  if (dim > 0) delete J;
}

template<class kind>
void Solver<kind>::set_default_values()
{
  max_its = 100;
  epsilon = 0.00001;
  dim = 0;
  nnonzero = 0;
  method = DIRECT;
}

namespace SYNARMOSMA {
  template<>
  double Solver<double>::a1(double t) const
  {
    // The first complement to the homotopy function,
    // s(t) = a1(t)*F(x) + a2(t)*x
    // We have the conditions on a1(t) (0 <= t <= 1) that 
    // it is continuous and
    // a1(0) = 0
    // a1(1) = 1
    double output = t;
    return output;
  }

  template<>
  std::complex<double> Solver<std::complex<double> >::a1(double t) const
  {
    // The first complement to the homotopy function,
    // s(t) = a1(t)*F(x) + a2(t)*x
    // We have the conditions on a1(t) (0 <= t <= 1) that 
    // it is continuous and
    // a1(0) = 0
    // a1(1) = 1
    std::complex<double> i(0.0,1.0);
    std::complex<double> output = std::sin(M_PI/2.0*t + 5.0*i*t*(t - 1.0));
    return output;
  }

  template<>
  double Solver<double>::a2(double t) const
  {
    // The second complement to the homotopy function,
    // s(t) = a1(t)*F(x) + a2(t)*x
    // We have the conditions on a2(t) (0 <= t <= 1) that 
    // it is continuous and
    // a2(0) = 1
    // a2(1) = 0
    double output = 1.0 - t;
    return output;
  }

  template<>
  std::complex<double> Solver<std::complex<double> >::a2(double t) const
  {
    // The second complement to the homotopy function,
    // s(t) = a1(t)*F(x) + a2(t)*x
    // We have the conditions on a2(t) (0 <= t <= 1) that 
    // it is continuous and
    // a2(0) = 1
    // a2(1) = 0
    std::complex<double> i(0.0,1.0);
    std::complex<double> output = (1.0 - t)*std::exp(10.0*i*t);
    return output;
  }
}

template<class kind>
void Solver<kind>::compute_jacobian(const std::vector<kind>& x)
{
  unsigned int i,j;
  static bool fcall = true;
  kind delta;
  std::set<unsigned int>::const_iterator it;
  std::vector<kind> w,y1,y2;

  J->clear(false);
  F(x,y1);
  w = x;
  if (fcall) {
    compute_dependencies();
    fcall = false;
  }

  for(i=0; i<dim; ++i) {
    for(it=dependencies[i].begin(); it!=dependencies[i].end(); ++it) {
      j = *it;
      w[j] += epsilon;
      F(w,y2);
      delta = y2[j] - y1[j];
      J->set(i,j,delta/epsilon);
      w[j] -= epsilon;
    }
  }
}

template<class kind>
void Solver<kind>::compute_dependencies()
{
  unsigned int i,j;
  kind delta;
  std::set<unsigned int> null;
  std::vector<kind> w,x,y;

  dependencies.clear();
  for(i=0; i<dim; ++i) {
    dependencies.push_back(null);
    delta = RND.drandom(-5.0,5.0);
    x.push_back(delta);
    y.push_back(kind(0.0));
    w.push_back(kind(0.0));
  }

  F(x,w);
  for(i=0; i<dim; ++i) {
    x[i] += kind(10.0);
    F(x,y);
    for(j=0; j<dim; ++j) {
      delta = y[j] - w[j];
      if (std::abs(delta) > epsilon) dependencies[j].insert(i);
    }
    x[i] -= kind(10.0);
  }
}

namespace SYNARMOSMA {
  template<>
  bool Solver<double>::direct_solver(std::vector<double>& x) const
  {
    int info,one = 1,n = dim;
    int pivots[dim];
    double A[dim*dim];

    J->convert(A);

    dgesv_(&n,&one,A,&n,pivots,&x[0],&n,&info);

    if (info != 0) {
      return false;
    }
    else {
      return true;
    }
  }

  template<>
  bool Solver<std::complex<double> >::direct_solver(std::vector<std::complex<double> >& x) const
  {
    int info,one = 1,n = dim;
    int pivots[dim];
    std::complex<double> A[dim*dim];

    J->convert(A);

    zgesv_(&n,&one,A,&n,pivots,&x[0],&n,&info);

    if (info != 0) {
      return false;
    }
    else {
      return true;
    }
  }
}

template<class kind>
bool Solver<kind>::linear_solver(const std::vector<kind>& x,const std::vector<kind>& b,std::vector<kind>& xnew) const
{
  bool success;
  if (method == DIRECT) {
    // Use the direct LAPACK-based linear solver
    xnew = b;
    success = direct_solver(xnew);
  }
  else {
    // Use the native Gauss-Seidel iterative solver in the Matrix class
    xnew = x;
    int n = J->gauss_seidel_solver(xnew,b,epsilon,100);
    success = (n > 0) ? true : false;
  }
  return success;
}

template<class kind>
bool Solver<kind>::forward_step()
{
  unsigned int i,its = 0;
  bool success = false;
  double q,norm,dt = 0.05;
  std::vector<kind> x,xnew,y,z;
  const kind w1 = a1(t);
  const kind w2 = a2(t);

  for(i=0; i<dim; ++i) {
    x.push_back(c_solution[i]);
    xnew.push_back(kind(0.0));
    y.push_back(kind(0.0));
  }

  do {
    compute_jacobian(x);
    // Now compute the rhs of the equation, b = J*x - F(x)
    J->multiply(x,z);
    for(i=0; i<dim; ++i) {
      y[i] = z[i] - y[i];
    }
    // And solve the linear system J*xnew = b
    success = linear_solver(x,y,xnew);
    if (!success) break;
    norm = 0.0;
    for(i=0; i<dim; ++i) {
      q = std::abs(xnew[i] - x[i]);
      norm += q*q;
      x[i] = xnew[i];  
    }
    // Check to see if the new solution differs appreciably from 
    // the old one
    if (std::sqrt(norm) < epsilon) break;
    F(x,y);
    norm = 0.0;
    for(i=0; i<dim; ++i) {
      q = std::abs(w1*y[i] + w2*x[i]);
      norm += q*q;
    }
    // Check to see if the new solution in fact solves the nonlinear
    // system
    if (norm < epsilon) {
      success = true;
      break;
    }
  } while(its <= max_its);
  if (success) {
    for(i=0; i<dim; ++i) {
      c_solution[i] = xnew[i];
    }
    t += dt;
  }
  return success;
}

template<class kind>
bool Solver<kind>::solve(std::vector<kind>& output)
{
  unsigned int i;
  bool success;

  for(i=0; i<dim; ++i) {
    c_solution.push_back(kind(0.0));
  }
  output = c_solution;

  do {
    success = forward_step();
    if (!success) break;
    if (std::abs(1.0 - t) < epsilon) break;
  } while(true);

  if (success) output = c_solution;
  
  return success;
}
