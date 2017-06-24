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
  homotopy = false;
  broyden = false;
  method = ITERATIVE; 
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
int Solver<kind>::forward_step()
{
  unsigned int i,j,its = 0;
  bool success,output = false;
  double q,old_norm,xnorm,fnorm,pfactor;
  kind p;
  std::vector<kind> x,xnew,f,fnew,b,fdiff,xdiff,z;
  Matrix<kind> secant(dim);

  x = c_solution;
  for(i=0; i<dim; ++i) {
    b.push_back(kind(0.0));
    f.push_back(kind(0.0));
    fnew.push_back(kind(0.0));
    xnew.push_back(kind(0.0));
    fdiff.push_back(kind(0.0));
    xdiff.push_back(kind(0.0));
  }
  homotopy_function(x,f);
  compute_jacobian(x);

  old_norm = 0.0;
  for(i=0; i<dim; ++i) {
    q = std::abs(f[i]);
    old_norm += q*q;
  }
  old_norm = std::sqrt(old_norm);
#ifdef VERBOSE
  std::cout << "Initial function norm is " << old_norm << std::endl;
#endif
  do {
    // Now compute the rhs of the equation, b = J*x - F(x)
    J->multiply(x,z);
    for(i=0; i<dim; ++i) {
      b[i] = z[i] - f[i];
    }
    // And solve the linear system J*xnew = b
    success = linear_solver(x,b,xnew);
    if (!success) break;
    xnorm = 0.0;
    for(i=0; i<dim; ++i) {
      xdiff[i] = xnew[i] - x[i];
      q = std::abs(xdiff[i]);
      xnorm += q*q;
    }    
    xnorm = std::sqrt(xnorm);
    // Check to see if the new solution differs appreciably from 
    // the old one
#ifdef VERBOSE
    std::cout << "Iterative difference is " << xnorm << " at " << its << std::endl;
#endif
    if (std::sqrt(xnorm) < epsilon) break;
    homotopy_function(xnew,fnew);
    fnorm = 0.0;
    for(i=0; i<dim; ++i) {
      fdiff[i] = fnew[i] - f[i];
      q = std::abs(fnew[i]);
      fnorm += q*q;
    }
    fnorm = std::sqrt(fnorm);
    // Check to see if the new solution in fact solves the nonlinear
    // system
#ifdef VERBOSE
    std::cout << "Function norm is " << fnorm << " at " << its << std::endl;
#endif
    if (fnorm < epsilon) {
      output = true;
      break;
    }
    // Diverging, time to exit...
    if (std::abs(fnorm/old_norm) > 10.0) break;
    old_norm = fnorm;
    its++;
    if (its > max_its) break;
    if (broyden) {
      // Compute the Broyden approximation to the Jacobian and use it:
      J->multiply(xdiff,z);
      for(i=0; i<dim; ++i) {
        z[i] = fdiff[i] - z[i];
      }
      pfactor = 1.0/(xnorm*xnorm);
      secant.clear(false);
      for(i=0; i<dim; ++i) {
        for(j=0; j<dim; ++j) {
          p = z[j]*xdiff[i];
          if (std::abs(p) > std::numeric_limits<double>::epsilon()) secant.set(i,j,pfactor*p);
        }
      }
      J->increment(secant);
    }
    else {
      compute_jacobian(xnew);
    }
    x = xnew;
    f = fnew;
  } while(true);

  if (output) {
    c_solution = xnew;
    return its;
  }
  return -1;
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
  int n;
  double dt = 0.01;
  bool success = true;

  if (output.size() == dim) {
    base_solution = output;
  }
  else {
    initialize_base_solution();
    output = base_solution;
  }

  c_solution = base_solution;

  if (homotopy) {
    t = 0.01;
    do {
      n = forward_step();
#ifdef VERBOSE
      std::cout << "Homotopy progress: " << t << "  " << n << std::endl;
#endif
      if (n > 0 && std::abs(1.0 - t) < epsilon) break;
      if (n == -1) {
        dt /= 2.0;
      }
      else {
        if (n < int(0.25*max_its)) {
          dt *= 2.0;
        }
        else if (n > int(0.75*max_its)) {
          dt /= 2.0;
        }
        t += dt;
      }
      if (dt < epsilon) {
#ifdef VERBOSE
        std::cout << "Homotopy method failed, exiting..." << std::endl;
#endif
        success = false;
        break;
      }
      if (t > 1.0) t = 1.0;
    } while(true);
  } 
  else {
    t = 1.0;
    n = forward_step();
    if (n == -1) success = false;
  }

  if (success) output = c_solution;
  
  return success;
}
