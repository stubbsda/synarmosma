#include "random.h"
#include "matrix.h"

#ifndef _solverh
#define _solverh

namespace SYNARMOSMA {
  template<class kind>
  /// A template class representing a Newton-Raphson solver of a system of nonlinear algebraic equations over a floating point base type.
  class Solver {
   protected:
    enum class Linear_Solver
    {
        direct,
        iterative
    };

    /// This property stores the maximum number of 
    /// Newton-Raphson iterations that will be carried 
    /// out before the Solver declares the situation 
    /// to be non-convergent. 
    unsigned int max_its = 100;
    /// The number of equations in the system, which is 
    /// assumed to also be the number of independent variables. 
    unsigned int dimension = 0;
    /// The epislon used to determine when the system is solved, 
    /// so that the norm of the error vector is less than this value.
    double epsilon = 0.00001;  
    double t = 1.0;
    /// A Boolean property that determines whether or not a homotopy-based 
    /// solution method is used.
    bool homotopy = false;
    bool broyden = false;
    /// This enumerated class determines whether the Solver will use a direct 
    /// or iterative linear solver; the former calls an LAPACK routine and is 
    /// suitable for dense matrices, while the latter makes use of the Gauss-Seidel 
    /// solver in the Matrix class and is appropriate for sparse matrices that 
    /// diagonally dominant.
    Linear_Solver method = Linear_Solver::iterative;
    /// This STL vector holds the current or working solution of the 
    /// system.
    std::vector<kind> current_solution;
    std::vector<kind> base_solution;
    /// The vector of independent variable dependencies per equation; an individual 
    /// set in the vector corresponds to the independent variables that exist in that 
    /// equation. This can simplify dramatically the calculation of the Jacobian 
    /// matrix since most systems that arise from applications are very sparse, with 
    /// each individual equation containing only a handful of independent variables. 
    std::vector<std::set<unsigned int> > dependencies;
    /// The Jacobian matrix of the system, i.e. the matrix of partial derivatives 
    /// of each equation (the rows) with respect to each independent variable (the 
    /// columns).
    Matrix<kind>* J;

    int forward_step();
    bool direct_solver(std::vector<kind>&) const;
    /// This method computes the Jacobian matrix for a particular set of values of the independent variables, which is the method's unique argument.
    void compute_jacobian(const std::vector<kind>&);
    /// This method calculates the content of the Solver::dependencies vector and returns the density, i.e. the sum of the cardinalities divided by the square of Solver::dimension.
    double compute_dependencies();
    bool linear_solver(const std::vector<kind>&,const std::vector<kind>&,std::vector<kind>&) const;
    kind a1(double) const;
    kind a2(double) const;
    void initialize_base_solution();
    virtual void F(const std::vector<kind>&,std::vector<kind>&) const = 0;
   public:
    Solver(unsigned int);
    Solver(unsigned int,double,unsigned int,bool,bool);
    ~Solver();
    bool solve(std::vector<kind>&);
    inline void use_iterative() {method = Linear_Solver::iterative;};
    inline void use_direct() {method = Linear_Solver::direct;};
    inline void use_homotopy() {homotopy = true;};
    inline void use_broyden() {broyden = true;};
  };
}
#endif
