#include "random.h"
#include "graph.h"

#ifndef _solverh
#define _solverh

namespace SYNARMOSMA {
  template<class kind>
  /// An abstract template class representing a Newton-Raphson solver of a system of nonlinear algebraic equations over a floating point base type.
  class Solver {
   protected:
    /// This enumerated class lists the two methods available for solving the linear system at each 
    /// stage of the Newton-Raphson solver: a direct method (using a LAPACK routine) or a native iterative 
    /// solver from the Matrix class.
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
    /// This property stores the number of equations in the system, which is 
    /// assumed to also be the number of independent variables. 
    unsigned int dimension = 0;
    /// This property is used to determine when the system is solved, 
    /// so that the norm of the error vector is less than this value.
    double error_threshold = 0.00001;
    /// This floating point property is only used when Solver::homotopy is 
    /// true, in which case it tracks the progress of the homotopy from 
    /// \f$t=0\f$ to \f$t=1\f$.  
    double t = 1.0;
    /// A Boolean property that determines whether or not a homotopy-based 
    /// solution method is used. If true, the system to be solved becomes 
    /// \f$F_H(t,x) = a_1(t) F(x) + a_2(t) (x - b)\f$ where \f$b\f$ is the 
    /// Solver::base_solution property and \f$F(x)=0\f$ the system of equations 
    /// we wish to solve. We choose the function \f$a_1(t)\f$ to be continuous 
    /// for \f$0\le t\le 1\f$ with \f$a_1(0) = 0\f$ and \f$a_1(1)=1\f$, while 
    /// \f$a_2(t)\f$ is continuous on the same interval and satifies \f$a_2(0)=1\f$ 
    /// and \f$a_2(1)=0\f$. We then seek to gradually advance from \f$t=0\f$, 
    /// where \f$H(0,x) = x - b = 0\f$ (so that \f$x=b\f$) to \f$t=1\f$, where 
    /// \f$F_H(1,x) = F(x) = 0\f$. 
    bool homotopy = false;
    /// This Boolean property controls whether or not the Broyden approximation 
    /// to the Jacobian matrix is used. In this approach, at the first step 
    /// we calculate the Jacobian \f$J_0\f$ by the standard method of finite 
    /// differences but for \f$n\ge 1\f$ we use the secant approximation to 
    /// obtain the following expression for the Jacobian, \f$J_n = J_{n-1} + 
    /// [\Delta F_n - J_{n-1}\Delta x_n]/|\Delta x_n|^2 \Delta x_n^T\f$, where 
    /// \f$\Delta F_n = F(x_n) - F(x_{n-1})\f$ and \f$\Delta x_n = x - x_{n-1}\f$.
    /// Using the Broyden method should dramatically increase the speed of the 
    /// Jacobian calculation at the expense of having to perform more Newton-Raphson 
    /// iterations to satisfy the convergence threshold.  
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
    /// When using the homotopy method to solve the system of equations, this vector 
    /// stores the random solution used at the initial homotopy step \f$t = 0\f$.
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

    /// This method carries out an iterative Newton-Raphson solution of the nonlinear system, returning the number of such iterations carried out. The method throws a runtime error if the convergence fails for any reason, which is normally caught by the solve() method. 
    int forward_step();
    /// This method attempts to solve the linear system \f$Jx = b\f$ where \f$b\f$ is the argument to the method and \f$J\f$ is the Solver::J matrix. The argument contains the solution \f$x\f$ if the method returns true, meaning the system was successfully solved using an LAPACK function. 
    bool direct_solver(std::vector<kind>&) const;
    /// This method computes the Jacobian matrix for a particular set of values of the independent variables, which is the method's unique argument.
    void compute_jacobian(const std::vector<kind>&);
    /// This method calculates the content of the Solver::dependencies vector and returns the density of Solver::J, i.e. the sum of the cardinalities divided by the square of Solver::dimension.
    double compute_dependencies();
    /// This method is a generic entry point for the solution of \f$Jx = b\f$, with \f$J\f$ the Solver::J matrix, using either the direct_solver() method or the Gauss-Seidel method of the Matrix class. The first argument is the initial guess for the solution of this linear system, the second argument is the source term vector \f$b\f$ and the final argument is the computed solution, if the method returns true. 
    bool linear_solver(const std::vector<kind>&,const std::vector<kind>&,std::vector<kind>&) const;
    /// This method is only of importance if Solver::homotopy is true, in which case it is the coefficient function for the first (i.e. \f$F(x)\f$) term of the homotopy. It must be continuous on \f$0\le t\le 1\f$ and satisfy the two conditions \f$a_1(0)=0\f$ and \f$a_1(1)=1\f$.  
    kind a1(double) const;
    /// This method is only of importance if Solver::homotopy is true, in which case it is the coefficient function for the second (i.e. \f$x-b\f$) term of the homotopy. It must be continuous on \f$0\le t\le 1\f$ and satisfy the two conditions \f$a_2(0)=1\f$ and \f$a_2(1)=0\f$.  
    kind a2(double) const;
    /// This method is used when Solver::homotopy is true to initialize the contents of Solver::base_solution, the solution of the system when \f$t=0\f$. The contents of the vector are uniform random variates on the interval \f$[\rho-L,\rho+L]\f$ where \f$\rho\f$ is the first argument and \f$2L\f$ is the second argument. 
    void initialize_base_solution(double,double);
    /// This is a pure virtual function which defines the system of equations \f$F(x) = 0\f$ that must be solved - the first argument contains the input \f$x\f$ while the second argument is the output.
    virtual void F(const std::vector<kind>&,std::vector<kind>&) const = 0;
   public:
    /// This constructor accepts as its unique argument the value of Solver::dimension, with all other properties retaining their default value.
    Solver(unsigned int);
    /// This constructor allows maximum flexibility: the arguments are the value of Solver::dimension, Solver::epsilon, Solver::max_its, Solver::homotopy and Solver::broyden.
    Solver(unsigned int,double,unsigned int,bool,bool);
    /// The destructor frees the memory associated with Solver::J if Solver::dimension is greater than zero.
    ~Solver();
    /// This method restores all of the instance's scalar properties to their default value, deletes the Solver::J property and clears the other extended properties.
    void clear();
    /// This method is the main entry point for solving the system of equations; the unique argument is the initial guess (which only makes sense if Solver::homotopy is false) which upon exit will contain the solution, if the return value is true.
    bool solve(std::vector<kind>&);
    /// This method computes the dependency graph among the system of equations - each equation corresponds to a vertex and if two equations have at least one independent variable in common, there is an edge connecting the corresponding vertices. The method returns the connectedness of this graph.
    bool compute_dependency_graph(Graph*) const;
    /// This is a public method that sets the Solver::method property to Linear_Solver::iterative.
    void use_iterative();
    /// This is a public method that sets the Solver::method property to Linear_Solver::direct.
    void use_direct();
    /// This is a public method that sets the Solver::homotopy property to true.
    void use_homotopy();
    /// This is a public method that sets the Solver::broyden property to true.
    void use_broyden();
    /// This method writes the instance properties to a binary disk file and returns the number of bytes written to the file.
    int serialize(std::ofstream&) const;
    /// This method calls the clear() method of the instance and then reads the properties from a binary disk file, returning the number of bytes read.
    int deserialize(std::ifstream&);
  };

  template<class kind>
  inline void Solver<kind>::use_iterative() 
  {
    method = Linear_Solver::iterative;
  }

  template<class kind>
  inline void Solver<kind>::use_direct() 
  {
    method = Linear_Solver::direct;
  }

  template<class kind>
  inline void Solver<kind>::use_homotopy() 
  {
    homotopy = true;
  }

  template<class kind>
  inline void Solver<kind>::use_broyden() 
  {
    broyden = true;
  }
}
#endif
