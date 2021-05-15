#include "multitime.h"

using namespace SYNARMOSMA;

template<class kind>
Multitime<kind>::Multitime()
{
  allocate();
}

template<class kind>
Multitime<kind>::Multitime(kind t)
{
  allocate();
  chronos[0].first = t;
  for(int i=1; i<Multitime<kind>::tdimension; ++i) {
    chronos[i].first = kind(0);
    chronos[i].second = false;
  }
}

template<class kind>
Multitime<kind>::Multitime(const Multitime<kind>& source)
{
  allocate();
  for(int i=0; i<Multitime<kind>::tdimension; ++i) {
    chronos[i] = source.chronos[i];
  }
}

template<class kind>
Multitime<kind>& Multitime<kind>::operator=(const Multitime<kind>& source)
{
  if (this == &source) return *this;

  allocate();
  for(int i=0; i<Multitime<kind>::tdimension; ++i) {
    chronos[i] = source.chronos[i];
  }

  return *this;
}

template<class kind>
Multitime<kind>& Multitime<kind>::operator -(const Multitime<kind>& source)
{
  for(int i=0; i<Multitime<kind>::tdimension; ++i) {
    chronos[i] = source.chronos[i];
    if (chronos[i].second) chronos[i].first = -source.chronos[i].first;
  }

  return *this;
}

template<class kind>
Multitime<kind>::~Multitime()
{
  delete[] chronos;
}

template<class kind>
void Multitime<kind>::clear()
{
  for(int i=0; i<Multitime<kind>::tdimension; ++i) {
    chronos[i].first = kind(0);
    chronos[i].second = true;
  }
}

template<class kind>
void Multitime<kind>::allocate()
{
  chronos = new std::pair<kind,bool>[Multitime<kind>::tdimension];
  for(int i=0; i<Multitime::tdimension; ++i) {
    chronos[i] = std::pair<kind,bool>(kind(0),true);
  }
}

namespace SYNARMOSMA {
  template<>
  /// This method is an instantiation of set() for the case of a discrete time vector, which is needed because in this case the method's initial argument first needs to be converted to a 64 bit integer through division by Multitime::time_quantum.
  void Multitime<INT64>::set(double alpha,int n)
  {
    if (n < 0 || n >= Multitime<INT64>::tdimension) throw std::invalid_argument("Illegal Multitime::set index!");

    chronos[n].first = INT64(alpha/time_quantum);
    chronos[n].second = true;
    for(int i=0; i<Multitime<INT64>::tdimension; ++i) {
      if (i == n) continue;
      chronos[i].second = false;
    }
  }

  template<>
  /// This method is an instantiation of set() for the case of a discrete time vector, which is needed because in this case the elements of the input vector first need to be converted to 64 bit integers through division by Multitime::time_quantum.
  void Multitime<INT64>::set(const std::vector<double>& alpha)
  {
    const int n = (signed) alpha.size();
    if (n > Multitime<INT64>::tdimension) throw std::invalid_argument("Illegal vector length in Multitime::set method!");

    int i;
    for(i=0; i<n; ++i) {
      chronos[i].first = INT64(alpha[i]/time_quantum);
      chronos[i].second = true;
    }
    for(i=n; i<Multitime<INT64>::tdimension; ++i) {
      chronos[i].first = 0;
      chronos[i].second = false;
    }
  }
}

template<class kind>
void Multitime<kind>::set(double alpha,int n)
{
  if (n < 0 || n >= Multitime<kind>::tdimension) throw std::invalid_argument("Illegal Multitime::set index!");

  chronos[n].first = alpha;
  chronos[n].second = true;
  for(int i=0; i<Multitime<kind>::tdimension; ++i) {
    if (i == n) continue;
    chronos[i].second = false;
  }
}

template<class kind>
void Multitime<kind>::set(const std::vector<double>& alpha)
{
  const int n = (signed) alpha.size();
  if (n > Multitime<kind>::tdimension) throw std::invalid_argument("Illegal vector length in Multitime::set method!");

  int i;
  for(i=0; i<n; ++i) {
    chronos[i].first = alpha[i];
    chronos[i].second = true;
  }
  for(i=n; i<Multitime<kind>::tdimension; ++i) {
    chronos[i].first = kind(0);
    chronos[i].second = false;
  }
}

template<class kind>
int Multitime<kind>::serialize(std::ofstream& s) const
{
  int count = 0;
  bool t;
  kind x;

  s.write((char*)(&Multitime<kind>::tdimension),sizeof(int)); count += sizeof(int);
  for(int i=0; i<Multitime<kind>::tdimension; ++i) {
    x = chronos[i].first;
    s.write((char*)(&x),sizeof(kind)); count += sizeof(kind);
    t = chronos[i].second;
    s.write((char*)(&t),sizeof(bool)); count += sizeof(bool);
  }
  return count;
}

template<class kind>
int Multitime<kind>::deserialize(std::ifstream& s)
{
  int n = Multitime<kind>::tdimension,count = 0;
  bool t;
  kind x;

  s.read((char*)(&Multitime::tdimension),sizeof(int)); count += sizeof(int);
  if (n != Multitime::tdimension) {
    delete[] chronos;
    chronos = new std::pair<kind,bool>[Multitime<kind>::tdimension];
  }

  for(int i=0; i<Multitime::tdimension; ++i) {
    s.read((char*)(&x),sizeof(kind)); count += sizeof(kind);
    s.read((char*)(&t),sizeof(bool)); count += sizeof(bool);
    chronos[i] = std::pair<double,bool>(x,t);
  }
  return count;
}

namespace SYNARMOSMA {
  template<>
  /// This method is an instantiation of norm() for the case of a discrete time vector, which is needed because in this case the temporal elements first need to be converted to doubles through multiplication by Multitime::time_quantum.
  double Multitime<INT64>::norm() const
  {
    double alpha,sum = 0.0;
  
    for(int i=0; i<Multitime<INT64>::tdimension; ++i) {
      if (!chronos[i].second) continue;
      alpha = time_quantum*double(chronos[i].first);
      sum += alpha*alpha;
    }
    return sum;
  }

  template<>
  /// This method is an instantiation of compactify() for the case of a discrete time vector, which is needed because in this case the temporal elements first need to be converted to doubles through multiplication by Multitime::time_quantum.
  double Multitime<INT64>::compactify(double R) const
  {
    if (!chronos[0].second) std::invalid_argument("The initial time dimension must be active in the Multitime::compactify method!");
    double alpha,sum = time_quantum*double(chronos[0].first);
    const double pfactor = 2.0/M_PI;

    for(int i=1; i<Multitime<INT64>::tdimension; ++i) {
      if (!chronos[i].second) continue;
      alpha = time_quantum*double(chronos[i].first);
      alpha = pfactor*std::atan(alpha);
      sum += alpha/R;
    }
    return sum;
  }

  template<>
  /// This method is an instantiation of extract() for the case of a discrete time vector, which is needed because in this case the temporal elements first need to be converted to doubles through multiplication by Multitime::time_quantum.
  void Multitime<INT64>::extract(std::vector<double>& tau) const
  {
    double alpha;

    tau.clear();

    for(int i=0; i<Multitime<INT64>::tdimension; ++i) {
      if (!chronos[i].second) continue;
      alpha = time_quantum*double(chronos[i].first);
      tau.push_back(alpha);
    }
  }
}

template<class kind>
double Multitime<kind>::norm() const
{
  double sum = 0.0;
  
  for(int i=0; i<Multitime<kind>::tdimension; ++i) {
    if (!chronos[i].second) continue;
    sum += chronos[i].first*chronos[i].first;
  }
  return sum;
}

template<class kind>
double Multitime<kind>::compactify(double R) const
{
  if (!chronos[0].second) std::invalid_argument("The initial time dimension must be active in the Multitime::compactify method!");
  double alpha,sum = chronos[0].first;
  const double pfactor = 2.0/M_PI;

  for(int i=1; i<Multitime<kind>::tdimension; ++i) {
    if (!chronos[i].second) continue;
    alpha = pfactor*std::atan(chronos[i].first);
    sum += alpha/R;
  }
  return sum;
}

template<class kind>
void Multitime<kind>::extract(std::vector<double>& tau) const
{
  tau.clear();

  for(int i=0; i<Multitime<kind>::tdimension; ++i) {
    if (!chronos[i].second) continue;
    tau.push_back(chronos[i].first);
  }
}

