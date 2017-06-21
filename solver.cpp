#include "solver.h"

using namespace SYNARMOSMA;

extern Random RND;

template<class kind>
Solver<kind>::Solver(int N)
{
  assert(N > 0);
  set_default_values();
  dim = N;
  J = new SYNARMOSMA::Matrix<kind>(dim);
}

template<class kind>
Solver<kind>::Solver(int N,double eps,int M,bool htype,bool approx)
{
  assert(N > 0 && M > 0 && eps > 0.0);
  set_default_values();
  dim = N;
  epsilon = eps;
  max_its = M;
  homotopy = htype;
  broyden = approx;
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
  homotopy = true;
  broyden = false;
  method = ITERATIVE; 
  //method = DIRECT;
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
void Solver<kind>::homotopy_function(const std::vector<kind>& x,std::vector<kind>& output) const
{
  unsigned int i;
  kind q;
  std::vector<kind> y;
  const kind w1 = a1(t);
  const kind w2 = a2(t);

  F(x,y);

  output.clear();
  for(i=0; i<dim; ++i) {
    q = w1*y[i] + w2*(x[i] - base_solution[i]);
    output.push_back(q);
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
  homotopy_function(x,y1);
  w = x;
  if (fcall) {
    compute_dependencies();
    fcall = false;
  }

  for(i=0; i<dim; ++i) {
    for(it=dependencies[i].begin(); it!=dependencies[i].end(); ++it) {
      j = *it;
      w[j] += epsilon;
      homotopy_function(w,y2);
      delta = y2[i] - y1[i];
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

  homotopy_function(x,w);
  for(i=0; i<dim; ++i) {
    x[i] += kind(10.0);
    homotopy_function(x,y);
    for(j=0; j<dim; ++j) {
      delta = y[j] - w[j];
      if (std::abs(delta) > epsilon) dependencies[j].insert(i);
    }
    x[i] -= kind(10.0);
  }
#ifdef VERBOSE
  for(i=0; i<dim; ++i) {
    std::cout << "For equation " << 1+i << " there are " << dependencies[i].size() << " independent variables" << std::endl;
  }
#endif
}

namespace SYNARMOSMA {
  template<>
  bool Solver<double>::direct_solver(std::vector<double>& x) const
  {
    int info,one = 1,n = dim;
    int pivots[dim];
    double A[dim*dim];
    bool output = false;

    J->convert(A,'c');

    dgesv_(&n,&one,A,&n,pivots,&x[0],&n,&info);

    if (info == 0) output = true;

    return output;
  }

  template<>
  bool Solver<std::complex<double> >::direct_solver(std::vector<std::complex<double> >& x) const
  {
    int info,one = 1,n = dim;
    int pivots[dim];
    std::complex<double> A[dim*dim];
    bool output = false;

    J->convert(A,'c');

    zgesv_(&n,&one,A,&n,pivots,&x[0],&n,&info);

    if (info == 0) output = true;

    return output;
  }
}

template<class kind>
bool Solver<kind>::linear_solver(const std::vector<kind>& x,const std::vector<kind>& b,std::vector<kind>& xnew) const
{
  bool success = false;
  if (method == DIRECT) {
    // Use the direct LAPACK-based linear solver
    xnew = b;
    success = direct_solver(xnew);
  }
  else {
    // Use the native Gauss-Seidel iterative solver in the Matrix class
    xnew = x;
    int n = J->gauss_seidel_solver(xnew,b,epsilon,100);
    if (n > 0) success = true;
  }
  return success;
}

template<class kind>
bool Solver<kind>::forward_step()
{
  unsigned int i,its = 0;
  bool success,output = false;
  double q,old_norm,norm,dt = 0.05;
  std::vector<kind> x,xnew,b,y,z;
  Matrix<kind> secant(dim);

  for(i=0; i<dim; ++i) {
    x.push_back(c_solution[i]);
    xnew.push_back(kind(0.0));
    y.push_back(kind(0.0));
    b.push_back(kind(0.0));
  }
  homotopy_function(x,y);
  old_norm = 0.0;
  for(i=0; i<dim; ++i) {
    q = std::abs(y[i]);
    old_norm += q*q;
  }
  old_norm = std::sqrt(old_norm);
  compute_jacobian(x);

  do {
    // Now compute the rhs of the equation, b = J*x - F(x)
    J->multiply(x,z);
    for(i=0; i<dim; ++i) {
      b[i] = z[i] - y[i];
    }
    // And solve the linear system J*xnew = b
    success = linear_solver(x,b,xnew);
    if (!success) break;
    norm = 0.0;
    for(i=0; i<dim; ++i) {
      q = std::abs(xnew[i] - x[i]);
      norm += q*q;
    }
    norm = std::sqrt(norm);
    // Check to see if the new solution differs appreciably from 
    // the old one
#ifdef VERBOSE
    std::cout << "Iterative difference is " << norm << " at " << its << std::endl;
#endif
    if (std::sqrt(norm) < epsilon) break;
    x = xnew;
    homotopy_function(x,y);
    norm = 0.0;
    for(i=0; i<dim; ++i) {
      q = std::abs(y[i]);
      norm += q*q;
    }
    norm = std::sqrt(norm);
    // Check to see if the new solution in fact solves the nonlinear
    // system
#ifdef VERBOSE
    std::cout << "Function norm is " << norm << " at " << its << std::endl;
#endif
    if (norm < epsilon) {
      output = true;
      break;
    }
    // Diverging, time to exit...
    if (std::abs(norm/old_norm) > 10.0) break;
    old_norm = norm;
    its++;
    if (its > max_its) break;
    if (broyden) {
      // Use the Broyden approximation to the Jacobian:
      // J_n = J_{n-1} + 1/norm(\Delta x)^2 * (\Delta F - J_{n-1}\Delta x)\outer \Delta x
      J->increment(secant);
    }
    else {
      compute_jacobian(x);
    }
  } while(true);
  if (output) {
    c_solution = xnew;
    t += dt;
  }
  return output;
}

namespace SYNARMOSMA {
  template<>
  void Solver<double>::initialize_base_solution()
  {
    unsigned int i;
    double alpha;

    base_solution.clear();

    for(i=0; i<dim; ++i) {
      alpha = -1.0 + 2.0*RND.drandom();
      base_solution.push_back(alpha);
    }
  }

  template<>
  void Solver<std::complex<double> >::initialize_base_solution()
  {
    unsigned int i;
    double alpha,beta;

    base_solution.clear();

    for(i=0; i<dim; ++i) {
      alpha = -1.0 + 2.0*RND.drandom();
      beta = -1.0 + 2.0*RND.drandom();
      base_solution.push_back(std::complex<double>(alpha,beta));
    }
  }
}

template<class kind>
bool Solver<kind>::solve(std::vector<kind>& output)
{
  bool success;

  initialize_base_solution();

  c_solution = base_solution;
  output = c_solution;

  if (homotopy) {
    t = 0.01;
    do {
      success = forward_step();
#ifdef VERBOSE
      std::cout << "Homotopy progress: " << t << "  " << success << std::endl;
#endif
      if (!success) break;
      if (std::abs(1.0 - t) < epsilon) break;
    } while(true);
  } 
  else {
    t = 1.0;
    c_solution[0] = 3.5;
    c_solution[1] = -0.5;
    c_solution[2] = 0.1;
    success = forward_step();
  }

  if (success) output = c_solution;
  
  return success;
}
