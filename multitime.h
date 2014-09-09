#include "global.h"

#ifndef _mtimeh
#define _mtimeh

class Multitime {
 private:
  double v1[tdimension],v2[tdimension];
  bool active[tdimension];

  void clear();
 public:
  Multitime();
  Multitime(double);
  Multitime(const Multitime&);
  ~Multitime();
  Multitime& operator =(const Multitime&);
  void initialize(double);
  double norm() const;
  void extract(std::vector<double>&) const;
  friend Multitime operator +(const Multitime&,const Multitime&);
  friend Multitime operator -(const Multitime&,const Multitime&);
  friend Multitime operator *(const Multitime&,const Multitime&);
  friend bool operator >(const Multitime&,const Multitime&);
  friend bool operator <(const Multitime&,const Multitime&);
  friend bool operator ==(const Multitime&,const Multitime&);
  friend bool operator !=(const Multitime&,const Multitime&);
  friend bool operator >=(const Multitime&,const Multitime&);
  friend bool operator <=(const Multitime&,const Multitime&);
  friend Multitime operator *(double alpha,const Multitime&);
};
#endif 
