/*
  Copyright 2014 Daniel Stubbs

  This file is part of Synarmosma.

  Synarmosma is free software: you can redistribute it and/or modify 
  it under the terms of the GNU General Public License as published by 
  the Free Software Foundation, either version 3 of the License, or 
  (at your option) any later version.

  Synarmosma is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with Synarmosma.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "multitime.h"

using namespace SYNARMOSMA;

const unsigned int Multitime::tdimension;

Multitime::Multitime()
{
  allocate();
}

Multitime::Multitime(double t)
{
  unsigned int i;
  allocate();
  chronos[0].first = t;
  chronos[0].second = true;
  for(i=1; i<Multitime::tdimension; ++i) {
    chronos[i].second = false;
  }
}

Multitime::Multitime(const Multitime& source)
{
  chronos = source.chronos;
}

Multitime& Multitime::operator=(const Multitime& source)
{
  if (this != &source) 
  chronos = source.chronos;
  return *this;
}

Multitime::~Multitime()
{
  chronos.clear();
}

void Multitime::clear()
{
  unsigned int i;
  for(i=0; i<Multitime::tdimension; ++i) {
    chronos[i].first = 0.0;
    chronos[i].second = true;
  }
}

void Multitime::allocate()
{
  unsigned int i;
  for(i=0; i<Multitime::tdimension; ++i) {
    chronos.push_back(std::pair<double,bool>(0.0,true));
  }
}

void Multitime::initialize(double alpha)
{
  unsigned int i;
  chronos[0].first = alpha;
  chronos[1].second = true;
  for(i=1; i<Multitime::tdimension; ++i) {
    chronos[i].second = false;
  }
}

void Multitime::serialize(std::ofstream& s) const
{
  unsigned int i,n = chronos.size();
  bool t;
  double x;

  s.write((char*)(&Multitime::tdimension),sizeof(int));
  s.write((char*)(&n),sizeof(int));
  for(i=0; i<n; ++i) {
    x = chronos[i].first;
    s.write((char*)(&x),sizeof(double));
    t = chronos[i].second;
    s.write((char*)(&t),sizeof(bool));
  }
}

void Multitime::deserialize(std::ifstream& s)
{
  unsigned int i,n;
  bool t;
  double x;

  chronos.clear();

  s.read((char*)(&Multitime::tdimension),sizeof(int));
  s.read((char*)(&n),sizeof(int));
  for(i=0; i<n; ++i) {
    s.read((char*)(&x),sizeof(double));
    s.read((char*)(&t),sizeof(bool));
    chronos.push_back(std::pair<double,bool>(x,t));
  }
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
    unsigned int i;
    Multitime output;

    for(i=0; i<Multitime::tdimension; ++i) {
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
    unsigned int i;
    Multitime output;

    for(i=0; i<Multitime::tdimension; ++i) {
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
    unsigned int i;
    Multitime output = tau;

    for(i=0; i<Multitime::tdimension; ++i) {
      if (output.chronos[i].second) output.chronos[i].first *= alpha;
    }
    return output;
  }

  Multitime operator *(const Multitime& t1,const Multitime& t2)
  {
    unsigned int i;
    Multitime output;

    for(i=0; i<Multitime::tdimension; ++i) {
      if (t1.chronos[i].second && t2.chronos[i].second) {
        output.chronos[i].first = t1.chronos[i].first * t2.chronos[i].first;
      }
      else {
        output.chronos[i].second = false;
      }
    }
    return output;
  }
}

