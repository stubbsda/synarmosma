#include "solver.h"

using namespace SYNARMOSMA;

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
}

template<class kind>
kind Solver<kind>::a1(double t) const
{
  // The first complement to the homotopy function,
  // s(t) = a1(t)*F(x) + a2(t)*x
  // We have the conditions on a1(t) (0 <= t <= 1) that 
  // it is continuous and
  // a1(0) = 0
  // a1(1) = 1
  kind output = t;
  return output;
}

template<class kind>
kind Solver<kind>::a2(double t) const
{
  // The second complement to the homotopy function,
  // s(t) = a1(t)*F(x) + a2(t)*x
  // We have the conditions on a2(t) (0 <= t <= 1) that 
  // it is continuous and
  // a2(0) = 1
  // a2(1) = 0
  kind output = 1.0 - t;
  return output;
}

template<class kind>
void Solver<kind>::compute_jacobian(const std::vector<kind>& x)
{
  unsigned int i,j;
  kind delta;
  std::vector<kind> w,y1,y2;

  J->clear(false);
  F(x,y1);
  w = x;
  for(i=0; i<dim; ++i) {
    w[i] += epsilon;
    F(w,y2);
    for(j=0; j<dim; ++j) {
      delta = y2[j] - y1[j];
      if (std::abs(delta) > epsilon) J->set(i,j,delta/epsilon);
    }
    w[i] -= epsilon;
  }
}

template<class kind>
void Solver<kind>::linear_solver(const std::vector<kind>& x,const std::vector<kind>& b,std::vector<kind>& xnew) const
{
#ifdef LAPACK

#elif PETSC

#else
  // Use the native solver
#endif
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
    linear_solver(x,y,xnew);
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
  do {
    success = forward_step();
    if (!success) break;
    if (std::abs(1.0 - t) < epsilon) break;
  } while(true);
  if (success) {
    output.clear();
    for(i=0; i<dim; ++i) {
      output.push_back(c_solution[i]);
    }
  }
  return success;
}
