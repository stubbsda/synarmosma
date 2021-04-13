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
  for(int i=0; i<Multitime<kind>::tdimension; ++i) {
    chronos[i] = source.chronos[i];
  }
}

template<class kind>
Multitime<kind>& Multitime<kind>::operator=(const Multitime<kind>& source)
{
  if (this == &source) return *this;
 
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

template<class kind>
void Multitime<kind>::set(kind alpha,int n)
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

template<class kind>
double Multitime<kind>::norm() const
{
  double sum = 0.0;
  
  for(int i=0; i<Multitime<kind>::tdimension; ++i) {
    if (!chronos[i].second) continue;
    sum += double(chronos[i].first)*double(chronos[i].first);
  }
  return sum;
}

template<class kind>
void Multitime<kind>::extract(std::vector<kind>& tau) const
{
  tau.clear();

  for(int i=0; i<Multitime<kind>::tdimension; ++i) {
    if (!chronos[i].second) continue;
    tau.push_back(chronos[i].first);
  }
}

