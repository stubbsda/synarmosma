#include "multitime.h"

Multitime::Multitime()
{
  for(int i=0; i<tdimension; ++i) {
    active[i] = false;
  }
}

Multitime::Multitime(double t)
{
  v1[0] = t; v2[0] = t;
  active[0] = true;
  for(int i=1; i<tdimension; ++i) {
    active[i] = false;
  }
}

Multitime::Multitime(const Multitime& source)
{
  for(int i=0; i<tdimension; ++i) {
    v1[i] = source.v1[i];
    v2[i] = source.v2[i];
    active[i] = source.active[i];
  }
}

Multitime& Multitime::operator=(const Multitime& source)
{
  if (this != &source) {
    for(int i=0; i<tdimension; ++i) {
      v1[i] = source.v1[i];
      v2[i] = source.v2[i];
      active[i] = source.active[i];
    }
  }
  return *this;
}

Multitime::~Multitime()
{

}

void Multitime::clear()
{
  for(int i=0; i<tdimension; ++i) {
    active[i] = false;
  }
}

void Multitime::initialize(double alpha)
{
  v1[0] = alpha;
  v2[0] = alpha;
  active[0] = true;
  for(int i=1; i<tdimension; ++i) {
    active[i] = false;
  }
}

double Multitime::norm() const
{
  double sum = 0.0;
  
  for(int i=0; i<tdimension; ++i) {
    if (!active[i]) continue;
    sum += std::abs(v2[i] - v1[i]) + 0.5*(std::abs(v1[i]) + std::abs(v2[i]));
  }
  return sum;
}

void Multitime::extract(std::vector<double>& t) const
{
  t.clear();

  for(int i=0; i<tdimension; ++i) {
    if (!active[i]) continue;
    t.push_back(std::abs(v2[i] - v1[i]) + 0.5*(std::abs(v1[i]) + std::abs(v2[i])));
  }
}

bool operator <(const Multitime& t1,const Multitime& t2)
{
  // A difficult problem - how can we compare two \emph{multi-dimensional} times? 
  return true;
}

Multitime operator +(const Multitime& t1,const Multitime& t2)
{
  int i;
  Multitime output;

  output.clear();
  for(i=0; i<tdimension; ++i) {
    if (!t1.active[i] && !t2.active[i]) continue;
    output.active[i] = true;
    if (t1.active[i] && t2.active[i]) {
      output.v1[i] = t1.v1[i] + t2.v1[i];
      output.v2[i] = t1.v2[i] + t2.v2[i];
    }
    else if (t1.active[i] && !t2.active[i]) {
      output.v1[i] = t1.v1[i]; output.v2[i] = t1.v2[i];
    }
    else if (!t1.active[i] && t2.active[i]) {
      output.v1[i] = t2.v1[i]; output.v2[i] = t2.v2[i];
    }
  }
  return output;
}

