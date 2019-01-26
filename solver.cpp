#include "solver.h"

using namespace SYNARMOSMA;

extern Random RND;

template<class kind>
Solver<kind>::Solver(unsigned int N)
{
  if (N == 0) throw std::invalid_argument("The Solver::dimension property must be greater than zero!");
  dimension = N;
  J = new SYNARMOSMA::Matrix<kind>(dimension);
}

template<class kind>
Solver<kind>::Solver(unsigned int N,double eps,unsigned int M,bool htype,bool approx)
{
  if (N == 0) throw std::invalid_argument("The Solver::dimension property must be greater than zero!");
  if (M == 0) throw std::invalid_argument("The maximum number of iterations in the Solver class must be greater than zero!");
  if (eps < std::numeric_limits<double>::epsilon()) throw std::invalid_argument("The convergence threshold in the Solver class must be greater than zero!");

  dimension = N;
  epsilon = eps;
  max_its = M;
  homotopy = htype;
  broyden = approx;
  J = new SYNARMOSMA::Matrix<kind>(dimension);
}

template<class kind>
Solver<kind>::~Solver()
{
  if (dimension > 0) delete J;
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
  Matrix<kind>* JF = new Matrix<kind>(dimension);

  F(x,y1);
  w = x;
  if (fcall) {
    compute_dependencies();
    fcall = false;
  }

  for(i=0; i<dimension; ++i) {
    for(it=dependencies[i].begin(); it!=dependencies[i].end(); ++it) {
      j = *it;
      w[j] += epsilon;
      F(w,y2);
      delta = y2[i] - y1[i];
      JF->set(i,j,delta/epsilon);
      w[j] -= epsilon;
    }
  }
  JF->homotopy_scaling(a1(t),a2(t),J);
  delete JF;
}

template<class kind>
double Solver<kind>::compute_dependencies()
{
  unsigned int i,j;
  kind delta;
  std::set<unsigned int> null;
  std::vector<kind> w,x,y;

  dependencies.clear();
  for(i=0; i<dimension; ++i) {
    dependencies.push_back(null);
    delta = RND.drandom(-5.0,5.0);
    x.push_back(delta);
    y.push_back(kind(0.0));
    w.push_back(kind(0.0));
  }

  F(x,w);
  for(i=0; i<dimension; ++i) {
    x[i] += kind(10.0);
    F(x,y);
    for(j=0; j<dimension; ++j) {
      delta = y[j] - w[j];
      if (std::abs(delta) > epsilon) dependencies[j].insert(i);
    }
    x[i] -= kind(10.0);
  }

  unsigned int sum = 0;
  for(i=0; i<dimension; ++i) {
    sum += dependencies[i].size();
#ifdef VERBOSE
    std::cout << "For equation " << 1+i << " there are " << dependencies[i].size() << " independent variables" << std::endl;
#endif
  }
  return double(sum)/double(dimension*dimension);
}

template<class kind>
bool Solver<kind>::compute_dependency_graph(Graph* G) const
{
  // This method computes the dependency graph for the variety's equations.
  // There is a vertex for each equation and if two equations share at least
  // one independent variable there is an edge connecting the two vertices.
  unsigned int i,j;
  std::vector<unsigned int> v;

  G->clear();
  for(i=0; i<dimension; ++i) {
    G->add_vertex();
  }

  for(i=0; i<dimension; ++i) {
    for(j=1+i; j<dimension; ++j) {
      set_intersection(dependencies[i].begin(),dependencies[i].end(),dependencies[j].begin(),dependencies[j].end(),std::back_inserter(v));
      if (!v.empty()) G->add_edge(i,j);
      v.clear();
    }
  }
  bool output = G->connected();
  return output;
}

namespace SYNARMOSMA {
  template<>
  bool Solver<double>::direct_solver(std::vector<double>& x) const
  {
    int info,one = 1,n = dimension;
    int pivots[dimension];
    double A[dimension*dimension];
    bool output = false;

    J->convert(A,'c');

    dgesv_(&n,&one,A,&n,pivots,&x[0],&n,&info);

    if (info == 0) output = true;

    return output;
  }

  template<>
  bool Solver<std::complex<double> >::direct_solver(std::vector<std::complex<double> >& x) const
  {
    int info,one = 1,n = dimension;
    int pivots[dimension];
    std::complex<double> A[dimension*dimension];
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
  bool success = true;
  if (method == Linear_Solver::direct) {
    // Use the direct LAPACK-based linear solver
    xnew = b;
    success = direct_solver(xnew);
  }
  else {
    // Use the native Gauss-Seidel iterative solver in the Matrix class
    xnew = x;
    try {
      J->gauss_seidel_solver(xnew,b,epsilon,100);
    }
    catch (std::runtime_error& e) {
      success = false;
    }
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
  Matrix<kind> secant(dimension);

  x = current_solution;
  for(i=0; i<dimension; ++i) {
    b.push_back(kind(0.0));
    f.push_back(kind(0.0));
    fnew.push_back(kind(0.0));
    xnew.push_back(kind(0.0));
    fdiff.push_back(kind(0.0));
    xdiff.push_back(kind(0.0));
  }
  F(x,f);
  for(i=0; i<dimension; ++i) {
    f[i] = a1(t)*f[i] + a2(t)*(x[i] - base_solution[i]); 
  }
  compute_jacobian(x);

  old_norm = 0.0;
  for(i=0; i<dimension; ++i) {
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
    for(i=0; i<dimension; ++i) {
      b[i] = z[i] - f[i];
    }
    // And solve the linear system J*xnew = b
    success = linear_solver(x,b,xnew);
    if (!success) break;
    xnorm = 0.0;
    for(i=0; i<dimension; ++i) {
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
    if (xnorm < epsilon) break;
    F(xnew,fnew);
    for(i=0; i<dimension; ++i) {
      fnew[i] = a1(t)*fnew[i] + a2(t)*(xnew[i] - base_solution[i]); 
    }
    fnorm = 0.0;
    for(i=0; i<dimension; ++i) {
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
    // Diverging, time to break out of the principal loop...
    if ((fnorm/old_norm) > 10.0) throw std::runtime_error("Newton-Raphson iterations are diverging!");
    old_norm = fnorm;
    its++;
    if (its > max_its) break;
    if (broyden) {
      // Compute the Broyden approximation to the Jacobian and use it:
      J->multiply(xdiff,z);
      for(i=0; i<dimension; ++i) {
        z[i] = fdiff[i] - z[i];
      }
      pfactor = 1.0/(xnorm*xnorm);
      secant.clear(false);
      for(i=0; i<dimension; ++i) {
        for(j=0; j<dimension; ++j) {
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

  if (!output) throw std::runtime_error("Newton-Raphson solver failed to converge!");

  current_solution = xnew;
  return its;
}

namespace SYNARMOSMA {
  template<>
  void Solver<double>::initialize_base_solution(double mean,double range)
  {
    unsigned int i;
    const double lbound = mean - 0.5*range;
    const double ubound = mean + 0.5*range;

    base_solution.clear();

    for(i=0; i<dimension; ++i) {
      base_solution.push_back(RND.drandom(lbound,ubound));
    }
  }

  template<>
  void Solver<std::complex<double> >::initialize_base_solution(double mean,double range)
  {
    unsigned int i;
    double alpha,beta;
    const double lbound = mean - 0.5*range;
    const double ubound = mean + 0.5*range;

    base_solution.clear();

    for(i=0; i<dimension; ++i) {
      alpha = RND.drandom(lbound,ubound);
      beta = RND.drandom(lbound,ubound);
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

  if (output.size() == dimension) {
    base_solution = output;
  }
  else {
    initialize_base_solution(0.0,2.0);
    output = base_solution;
  }

  current_solution = base_solution;

  if (homotopy) {
    bool converged = false;
    t = 0.01;
    do {
      try {
        n = forward_step();
      }
      catch (std::runtime_error& e) {
#ifdef VERBOSE
        std::cout << "Homotopy solver failed at t = " << t << ", halving step size..." << std::endl;
#endif
        dt /= 2.0;
        t -= dt;
      }
      if (n > 0) {
        if (std::abs(1.0 - t) < epsilon) {
#ifdef VERBOSE
          std::cout << "Homotopy solver has converged, exiting..." << std::endl;
#endif
          converged = true;
        }
        else {
#ifdef VERBOSE
          std::cout << "Homotopy solver succeeded at t = " << t << " after " << n << " iterations." << std::endl;
#endif
          if (n < int(0.25*max_its)) {
            dt *= 2.0;
          }
          else if (n > int(0.75*max_its)) {
            dt /= 2.0;
          }
          t += dt;
        }
      }
      if (dt < epsilon) {
#ifdef VERBOSE
        std::cout << "Homotopy method has failed, exiting..." << std::endl;
#endif
        success = false;
        break;
      }
      if (t > 1.0) t = 1.0;
    } while(!converged);
  } 
  else {
    try {
      n = forward_step();
    }
    catch (std::runtime_error& e) {
      success = false;
    }
  }

  if (success) output = current_solution;
  
  return success;
}
