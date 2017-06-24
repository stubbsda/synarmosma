#include "matrix.h"

#ifndef _solverh
#define _solverh

namespace SYNARMOSMA {
  enum LINEAR_SOLVER_TYPE
  {
    DIRECT,
    ITERATIVE
  };

  template<class kind>
  class Solver {
   protected:
    unsigned int max_its;
    unsigned int dim;
    double epsilon;  
    double t;
    bool homotopy;
    bool broyden;
    LINEAR_SOLVER_TYPE method;
    std::vector<kind> c_solution,base_solution;
    std::vector<std::set<unsigned int> > dependencies;
    Matrix<kind>* J;

    int forward_step();
    bool direct_solver(std::vector<kind>&) const;
    void set_default_values();
    void compute_jacobian(const std::vector<kind>&);
    void compute_dependencies();
    bool linear_solver(const std::vector<kind>&,const std::vector<kind>&,std::vector<kind>&) const;
    kind a1(double) const;
    kind a2(double) const;
    void initialize_base_solution();
    void homotopy_function(const std::vector<kind>&,std::vector<kind>&) const;
    virtual void F(const std::vector<kind>&,std::vector<kind>&) const = 0;
   public:
    Solver(int);
    Solver(int,double,int,bool,bool);
    ~Solver();
    bool solve(std::vector<kind>&);
    inline void use_iterative() {method = ITERATIVE;};
    inline void use_direct() {method = DIRECT;};
    inline void use_homotopy() {homotopy = true;};
    inline void use_broyden() {broyden = true;};
  };
}
#endif
