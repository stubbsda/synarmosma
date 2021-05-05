#include "global.h"

#ifndef _mtimeh
#define _mtimeh

namespace SYNARMOSMA {
  template<class kind>
  class Multitime;

  template<class kind>
  std::ostream& operator <<(std::ostream&,const Multitime<kind>&);

  template<class kind>
  bool operator <=(const Multitime<kind>&,const Multitime<kind>&);

  template<class kind>
  bool operator <(const Multitime<kind>&,const Multitime<kind>&);

  template<class kind>
  bool operator >=(const Multitime<kind>&,const Multitime<kind>&);

  template<class kind>
  bool operator >(const Multitime<kind>&,const Multitime<kind>&);

  template<class kind>
  bool operator ==(const Multitime<kind>&,const Multitime<kind>&);

  template<class kind>
  bool operator !=(const Multitime<kind>&,const Multitime<kind>&);

  template<class kind>
  Multitime<kind> operator +(const Multitime<kind>&,const Multitime<kind>&);

  template<class kind>
  Multitime<kind> operator -(const Multitime<kind>&,const Multitime<kind>&);

  template<class kind>
  Multitime<kind> operator *(kind,const Multitime<kind>&);

  template<class kind>
  Multitime<kind> operator *(const Multitime<kind>&,const Multitime<kind>&);

  /// A class representing the concept of a dynamic multi-dimensional time.
  template<class kind>
  class Multitime {
   protected:
    /// This array represents the value of the time in each of its dimensions,
    /// along with a Boolean that indicates whether this temporal dimension is
    /// currently active.
    std::pair<kind,bool>* chronos;
    /// This property sets the length of the vector chronos and is so the maximum
    /// number of possible temporal dimensions; currently set to the mildly unconventional
    /// value of two.
    static const int tdimension = 2;
    /// This constant represents the smallest possible spatial separation between
    /// two moments of time; it is only meaningful when this template class is instantiated
    /// with a discrete base type.
    static const double time_quantum;

    /// This method sets each element of the Multitime::chronos array to (0.0,true), i.e. the same state as after the allocate() method.
    void clear();
    /// This method allocates the memory for the Multitime::chronos array and initializes each element to (0.0,true).
    void allocate();
   public:
    /// The default constructor, which just calls the allocate() method and exits.
    Multitime();
    /// This constructor calls the allocate method and then sets the first element to the argument and deactivates the other dimensions.
    Multitime(kind);
    /// The copy constructor just copies over the contents of the Multitime::chronos array.
    Multitime(const Multitime&);
    /// Destructor which frees the memory used by the Multitime::chronos array.
    ~Multitime();
    /// The assignment operator which copies over the contents of the Multitime::chronos array.
    Multitime& operator =(const Multitime&);
    /// This overloaded unary negation operator multiplies every element of the Multitime::chronos array by -1.
    Multitime& operator -(const Multitime&);
    /// This method writes the value of the Multitime::tdimension property and the contents of the Multitime::chronos array to a binary disk file, returning the number of bytes written.
    int serialize(std::ofstream&) const;
    /// This method first frees the memory of the Multitime::chronos array, reads in the Multitime::tdimension value, allocates the Multitime::chronos array and then reads in the contents of this array from a binary disk file and returns the number of bytes read.
    int deserialize(std::ifstream&);
    /// This method returns the value of the Multitime::tdimension property.
    static int get_dimension();
    /// This method returns the value of the sum of squares of the time values whose dimension is active.
    double norm() const;
    /// This method sets the time value of the element indicated by the second argument to the first argument and activates this element while deactivating all others.
    void set(kind,int = 0);
    /// This method extracts the time coordinate value for each active time dimension and puts it in a vector that will be the method's output.
    void extract(std::vector<kind>&) const;
    /// This operator assigns as output for each element the sum of the coordinate values if both argument dimensions are active, otherwise one or the other coordinative value if only one dimension is active.
    friend Multitime<kind> operator + <>(const Multitime<kind>&,const Multitime<kind>&);
    /// This operator assigns as output for each element the product of the coordinate values if both argument dimensions are active, otherwise it is zero and the dimension is inactive.
    friend Multitime<kind> operator * <>(const Multitime<kind>&,const Multitime<kind>&);
    /// This operator multiplies the coordinate of each active temporal dimension by the first argument.
    friend Multitime<kind> operator * <>(kind,const Multitime<kind>&);
    /// A Boolean operator that compares the norm of the two Multitime instances to determine their relative order.
    friend bool operator > <>(const Multitime<kind>&,const Multitime<kind>&);
    /// A Boolean operator that compares the norm of the two Multitime instances to determine their relative order.
    friend bool operator < <>(const Multitime<kind>&,const Multitime<kind>&);
    /// A Boolean operator that compares the norm of the two Multitime instances to determine their relative order.
    friend bool operator == <>(const Multitime<kind>&,const Multitime<kind>&);
    /// A Boolean operator that compares the norm of the two Multitime instances to determine their relative order.
    friend bool operator != <>(const Multitime<kind>&,const Multitime<kind>&);
    /// A Boolean operator that compares the norm of the two Multitime instances to determine their relative order.
    friend bool operator >= <>(const Multitime<kind>&,const Multitime<kind>&);
    /// A Boolean operator that compares the norm of the two Multitime instances to determine their relative order.
    friend bool operator <= <>(const Multitime<kind>&,const Multitime<kind>&);
    /// This overloading of the ostream operator writes the Multitime instance to the screen in a "pretty print" format.
    friend std::ostream& operator << <>(std::ostream&,const Multitime<kind>&);
  };

  template<class kind>
  inline int Multitime<kind>::get_dimension()
  {
    return tdimension;
  }

  template<class kind>
  std::ostream& operator <<(std::ostream& s,const Multitime<kind>& source)
  {
    int i,n = 0;

    for(i=0; i<Multitime<kind>::tdimension; ++i) {
      if (!source.chronos[i].second) continue;
      n++;
    }
    if (n > 0) { 
      int l = 0;

      s << "[";
      for(i=0; i<Multitime<kind>::tdimension; ++i) {
        if (!source.chronos[i].second) continue;
        l++;
        s << source.chronos[i].first;
        if (l < n) s << ",";
      }
      s << "]";
    }
    return s;
  }

  template<class kind>
  bool operator <(const Multitime<kind>& t1,const Multitime<kind>& t2)
  {
    // A difficult problem - how can we compare two \emph{multi-dimensional} times?
    return (t1.norm() < t2.norm());
  }

  template<class kind>
  bool operator >(const Multitime<kind>& t1,const Multitime<kind>& t2)
  {
    return (t1.norm() > t2.norm());
  }

  template<class kind>
  bool operator ==(const Multitime<kind>& t1,const Multitime<kind>& t2)
  {
    double alpha = std::abs(t1.norm() - t2.norm());
    if (alpha < std::numeric_limits<double>::epsilon()) return true;
    return false;
  }

  template<class kind>
  bool operator !=(const Multitime<kind>& t1,const Multitime<kind>& t2) 
  {
    return !(t1 == t2);
  }

  template<class kind>
  bool operator >=(const Multitime<kind>& t1,const Multitime<kind>& t2)
  {
    if ((t1 > t2) || (t1 == t2)) return true;
    return false;
  }

  template<class kind>
  bool operator <=(const Multitime<kind>& t1,const Multitime<kind>& t2)
  {
    if ((t1 < t2) || (t1 == t2)) return true;
    return false;
  }

  template<class kind>
  Multitime<kind> operator +(const Multitime<kind>& t1,const Multitime<kind>& t2)
  {
    Multitime<kind> output;

    for(int i=0; i<Multitime<kind>::tdimension; ++i) {
      output.chronos[i].second = false;
      if (!t1.chronos[i].second && !t2.chronos[i].second) continue;
      output.chronos[i].second = true;
      if (t1.chronos[i].second && t2.chronos[i].second) {
        output.chronos[i].first = t1.chronos[i].first + t2.chronos[i].first;
      }
      else if (t1.chronos[i].second && !t2.chronos[i].second) {
        output.chronos[i].first = t1.chronos[i].first;
      }
      else if (!t1.chronos[i].second && t2.chronos[i].second) {
        output.chronos[i].first = t2.chronos[i].first;
      }
    }
    return output;
  }

  template<class kind>
  Multitime<kind> operator -(const Multitime<kind>& t1,const Multitime<kind>& t2)
  {
    Multitime<kind> output;

    for(int i=0; i<Multitime<kind>::tdimension; ++i) {
      output.chronos[i].second = false;
      if (!t1.chronos[i].second && !t2.chronos[i].second) continue;
      output.chronos[i].second = true;
      if (t1.chronos[i].second && t2.chronos[i].second) {
        output.chronos[i].first = t1.chronos[i].first - t2.chronos[i].first;
      }
      else if (t1.chronos[i].second && !t2.chronos[i].second) {
        output.chronos[i].first = t1.chronos[i].first;
      }
      else if (!t1.chronos[i].second && t2.chronos[i].second) {
        output.chronos[i].first = -t2.chronos[i].first;
      }
    }
    return output;
  }

  template<class kind>
  Multitime<kind> operator *(kind alpha,const Multitime<kind>& tau)
  {
    Multitime<kind> output = tau;

    for(int i=0; i<Multitime<kind>::tdimension; ++i) {
      if (output.chronos[i].second) output.chronos[i].first *= alpha;
    }
    return output;
  }

  template<class kind>
  Multitime<kind> operator *(const Multitime<kind>& t1,const Multitime<kind>& t2)
  {
    Multitime<kind> output;

    for(int i=0; i<Multitime<kind>::tdimension; ++i) {
      if (t1.chronos[i].second && t2.chronos[i].second) {
        output.chronos[i].second = true;
        output.chronos[i].first = t1.chronos[i].first * t2.chronos[i].first;
      }
      else {
        output.chronos[i].second = false;
      }
    }
    return output;
  }
}
#endif 
