#include "group.h"

#ifndef _homotopyh
#define _homotopyh

class Homotopy {
 private:
  std::vector<Group> sequence;
  double fitness;

  void compute_fitness();
 public:
  Homotopy();
  Homotopy(int);
  Homotopy(const Homotopy&);
  ~Homotopy();
  void mutate();
  double get_fitness() const;
  friend Homotopy operator +(const Homotopy&,const Homotopy&);
  Homotopy& operator =(const Homotopy&);  
  void write();
};
#endif
