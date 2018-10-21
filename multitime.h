#include "global.h"

#ifndef _mtimeh
#define _mtimeh

namespace SYNARMOSMA {
  /// A class representing the concept of a dynamic multi-dimensional time.
  class Multitime {
   protected:
    /// This vector represents the value of the time in each of its dimensions, 
    /// along with a Boolean that indicates whether this temporal dimension is 
    /// currently active.
    std::vector<std::pair<double,bool> > chronos; 
    /// This property sets the length of the vector chronos and is so the maximum 
    /// number of possible temporal dimensions; currently set to the conventional 
    /// value of unity.
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
