#include "global.h"

#ifndef _mtimeh
#define _mtimeh

namespace SYNARMOSMA {
  class Multitime {
   protected:
    std::vector<std::pair<double,bool> > chronos; 
    // Not too Aristotelian, perhaps?
    static const unsigned int tdimension = 1;

    void clear();
    void allocate();
    void initialize(double);
   public:
    Multitime();
    Multitime(double);
    Multitime(const Multitime&);
    ~Multitime();
    Multitime& operator =(const Multitime&);
    int serialize(std::ofstream&) const;
    int deserialize(std::ifstream&);
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
    friend class Event;
  };
}
#endif 
