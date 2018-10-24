#include "multitime.h"

using namespace SYNARMOSMA;

const int Multitime::tdimension;

Multitime::Multitime()
{
  allocate();
}

Multitime::Multitime(double t)
{
  allocate();
  chronos[0].first = t;
  for(int i=1; i<Multitime::tdimension; ++i) {
    chronos[i].second = false;
  }
}

Multitime::Multitime(const Multitime& source)
{
  for(int i=0; i<Multitime::tdimension; ++i) {
    chronos[i] = source.chronos[i];
  }
}

Multitime& Multitime::operator=(const Multitime& source)
{
  if (this == &source) return *this;
 
  for(int i=0; i<Multitime::tdimension; ++i) {
    chronos[i] = source.chronos[i];
  }

  return *this;
}

Multitime::~Multitime()
{
  delete[] chronos;
}

void Multitime::clear()
{
  for(int i=0; i<Multitime::tdimension; ++i) {
    chronos[i].first = 0.0;
    chronos[i].second = true;
  }
}

void Multitime::allocate()
{
  chronos = new std::pair<double,bool>[Multitime::tdimension];
  for(int i=0; i<Multitime::tdimension; ++i) {
    chronos[i] = std::pair<double,bool>(0.0,true);
  }
}

void Multitime::set(double alpha,int n)
{
  if (n < 0 || n >= Multitime::tdimension) throw std::invalid_argument("Illegal Multitime::set index!");

  chronos[n].first = alpha;
  chronos[n].second = true;
  for(int i=0; i<Multitime::tdimension; ++i) {
    if (i == n) continue;
    chronos[i].second = false;
  }
}

int Multitime::serialize(std::ofstream& s) const
{
  int count = 0;
  bool t;
  double x;

  s.write((char*)(&Multitime::tdimension),sizeof(int)); count += sizeof(int);
  for(int i=0; i<Multitime::tdimension; ++i) {
    x = chronos[i].first;
    s.write((char*)(&x),sizeof(double)); count += sizeof(double);
    t = chronos[i].second;
    s.write((char*)(&t),sizeof(bool)); count += sizeof(bool);
  }
  return count;
}

int Multitime::deserialize(std::ifstream& s)
{
  int n = Multitime::tdimension,count = 0;
  bool t;
  double x;

  s.read((char*)(&Multitime::tdimension),sizeof(int)); count += sizeof(int);
  if (n != Multitime::tdimension) {
    delete[] chronos;
    chronos = new std::pair<double,bool>[Multitime::tdimension];
  }

  for(int i=0; i<Multitime::tdimension; ++i) {
    s.read((char*)(&x),sizeof(double)); count += sizeof(double);
    s.read((char*)(&t),sizeof(bool)); count += sizeof(bool);
    chronos[i] = std::pair<double,bool>(x,t);
  }
  return count;
}

double Multitime::norm() const
{
  double sum = 0.0;
  
  for(int i=0; i<Multitime::tdimension; ++i) {
    if (!chronos[i].second) continue;
    sum += chronos[i].first*chronos[i].first;
  }
  return sum;
}

void Multitime::extract(std::vector<double>& tau) const
{
  tau.clear();

  for(int i=0; i<Multitime::tdimension; ++i) {
    if (!chronos[i].second) continue;
    tau.push_back(chronos[i].first);
  }
}

namespace SYNARMOSMA {
  bool operator <(const Multitime& t1,const Multitime& t2)
  {
    // A difficult problem - how can we compare two \emph{multi-dimensional} times?
    return (t1.norm() < t2.norm());
  }

  bool operator >(const Multitime& t1,const Multitime& t2)
  {
    return (t1.norm() > t2.norm());
  }

  bool operator ==(const Multitime& t1,const Multitime& t2)
  {
    double alpha = std::abs(t1.norm() - t2.norm());
    if (alpha < std::numeric_limits<double>::epsilon()) return true;
    return false;
  }

  bool operator !=(const Multitime& t1,const Multitime& t2) 
  {
    return !(t1 == t2);
  }

  bool operator >=(const Multitime& t1,const Multitime& t2)
  {
    if ((t1 > t2) || (t1 == t2)) return true;
    return false;
  }

  bool operator <=(const Multitime& t1,const Multitime& t2)
  {
    if ((t1 < t2) || (t1 == t2)) return true;
    return false;
  }

  Multitime operator +(const Multitime& t1,const Multitime& t2)
  {
    Multitime output;

    for(int i=0; i<Multitime::tdimension; ++i) {
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

  Multitime operator -(const Multitime& t1,const Multitime& t2)
  {
    Multitime output;

    for(int i=0; i<Multitime::tdimension; ++i) {
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

  Multitime operator *(double alpha,const Multitime& tau)
  {
    Multitime output = tau;

    for(int i=0; i<Multitime::tdimension; ++i) {
      if (output.chronos[i].second) output.chronos[i].first *= alpha;
    }
    return output;
  }

  Multitime operator *(const Multitime& t1,const Multitime& t2)
  {
    Multitime output;

    for(int i=0; i<Multitime::tdimension; ++i) {
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

