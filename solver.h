#include "matrix.h"

#ifndef _solverh
#define _solverh

namespace SYNARMOSMA {
  template<class kind>
  class Solver {
   protected:
    unsigned int max_its;
    unsigned int dim;
    unsigned int nnonzero;
    double epsilon;  
    double t;
    std::vector<kind> c_solution;
    std::vector<std::set<unsigned int> > dependencies;
    Matrix<kind>* J;

    bool forward_step();
    bool direct_solver(std::vector<kind>&) const;
    void set_default_values();
    void compute_jacobian(const std::vector<kind>&);
    void compute_dependencies();
    bool linear_solver(const std::vector<kind>&,const std::vector<kind>&,std::vector<kind>&) const;
    kind a1(double) const;
    kind a2(double) const;
    virtual void F(const std::vector<kind>&,std::vector<kind>&) const = 0;
   public:
    Solver(int);
    ~Solver();
    bool solve(std::vector<kind>&);
  };
}
#endif
