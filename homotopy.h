#include "group.h"
#include "nexus.h"

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
  std::string write() const;
  void mutate();
  void clear();
  void compute(const Nexus*);
  void serialize(std::ofstream&) const;
  void deserialize(std::ifstream&);
  inline double get_fitness() const {return fitness;};
  friend Homotopy operator +(const Homotopy&,const Homotopy&);
  Homotopy& operator =(const Homotopy&);  
};
#endif
