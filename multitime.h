#include "global.h"

#ifndef _mtimeh
#define _mtimeh

namespace SYNARMOSMA {
  /// A class representing the concept of a dynamic multi-dimensional time.
  class Multitime {
   protected:
    /// This array represents the value of the time in each of its dimensions, 
    /// along with a Boolean that indicates whether this temporal dimension is 
    /// currently active.
    std::pair<double,bool>* chronos; 
    /// This property sets the length of the vector chronos and is so the maximum 
    /// number of possible temporal dimensions; currently set to the conventional 
    /// value of unity.
    static const int tdimension = 1;

    /// This method sets each element of the chronos array to (0.0,true), i.e. the same state as after the allocate method.
    void clear();
    /// This method allocates the memory for the chronos array and initializes each element to (0.0,true). 
    void allocate();
   public:
    /// The default constructor, which just calls the allocate method and exits.
    Multitime();
    /// This constructor calls the allocate method and then sets the first element to the argument and deactivates the other dimensions.
    Multitime(double);
    /// The copy constructor just copies over the contents of the chronos array.
    Multitime(const Multitime&);
    /// Destructor which frees the memory used by the chronos array.
    ~Multitime();
    /// The assignment operator which copies over the contents of the chronos array.
    Multitime& operator =(const Multitime&);
    /// This method writes the value of the tdimension property and the contents of the chronos array to a binary disk file, returning the number of bytes written.
    int serialize(std::ofstream&) const;
    /// This method first frees the memory of the chronos array, reads in the tdimension value, allocates the chronos array and then reads in the contents of this array from a binary disk file and returns the number of bytes read.
    int deserialize(std::ifstream&);
    /// This method returns the value of the sum of squares of the time values whose dimension is active.
    double norm() const;
    /// This method sets the time value of the element indicated by the second argument to the first argument and activates this element while deactivating all others.
    void set(double,int = 0);
    /// This method extracts the time coordinate value for each active time dimension and puts it in a vector that will be the method's output.
    void extract(std::vector<double>&) const;
    /// This operator assigns as output for each element the sum of the coordinate values if both argument dimensions are active, otherwise one or the other coordinative value if only one dimension is active.
    friend Multitime operator +(const Multitime&,const Multitime&);
    /// This operator assigns as output for each element the difference of the coordinate values if both argument dimensions are active, otherwise one or the other coordinative value if only one dimension is active.
    friend Multitime operator -(const Multitime&,const Multitime&);
    /// This operator assigns as output for each element the product of the coordinative values if both argument dimensions are active, otherwise it is zero and the dimension is inactive.
    friend Multitime operator *(const Multitime&,const Multitime&);
    /// A Boolean operator that compares the norm of the two Multitime instances to determine their relative order.
    friend bool operator >(const Multitime&,const Multitime&);
    /// A Boolean operator that compares the norm of the two Multitime instances to determine their relative order.
    friend bool operator <(const Multitime&,const Multitime&);
    /// A Boolean operator that compares the norm of the two Multitime instances to determine their relative order.
    friend bool operator ==(const Multitime&,const Multitime&);
    /// A Boolean operator that compares the norm of the two Multitime instances to determine their relative order.
    friend bool operator !=(const Multitime&,const Multitime&);
    /// A Boolean operator that compares the norm of the two Multitime instances to determine their relative order.
    friend bool operator >=(const Multitime&,const Multitime&);
    /// A Boolean operator that compares the norm of the two Multitime instances to determine their relative order.
    friend bool operator <=(const Multitime&,const Multitime&);
    /// This operator multiplies the coordinate of each active temporal dimension by the first argument.
    friend Multitime operator *(double alpha,const Multitime&);
  };
}
#endif 
