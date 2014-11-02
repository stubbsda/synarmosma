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

const int Multitime::tdimension;

Multitime::Multitime()
{
  allocate();
  for(int i=0; i<Multitime::tdimension; ++i) {
    active[i] = false;
  }
}

Multitime::Multitime(double t)
{
  allocate();
  v1[0] = t; v2[0] = t;
  active[0] = true;
  for(int i=1; i<Multitime::tdimension; ++i) {
    active[i] = false;
  }
}

Multitime::Multitime(const Multitime& source)
{
  allocate();
  for(int i=0; i<Multitime::tdimension; ++i) {
    v1[i] = source.v1[i];
    v2[i] = source.v2[i];
    active[i] = source.active[i];
  }
}

Multitime& Multitime::operator=(const Multitime& source)
{
  if (this != &source) {
    for(int i=0; i<Multitime::tdimension; ++i) {
      v1[i] = source.v1[i];
      v2[i] = source.v2[i];
      active[i] = source.active[i];
    }
  }
  return *this;
}

Multitime::~Multitime()
{
  delete[] v1;
  delete[] v2;
  delete[] active;
}

void Multitime::clear()
{
  for(int i=0; i<Multitime::tdimension; ++i) {
    active[i] = false;
  }
}

void Multitime::allocate()
{
  v1 = new double[Multitime::tdimension];
  v2 = new double[Multitime::tdimension];
  active = new bool[Multitime::tdimension];
}

void Multitime::initialize(double alpha)
{
  v1[0] = alpha;
  v2[0] = alpha;
  active[0] = true;
  for(int i=1; i<Multitime::tdimension; ++i) {
    active[i] = false;
  }
}

double Multitime::norm() const
{
  double sum = 0.0;
  
  for(int i=0; i<Multitime::tdimension; ++i) {
    if (!active[i]) continue;
    sum += std::abs(v2[i] - v1[i]) + 0.5*(std::abs(v1[i]) + std::abs(v2[i]));
  }
  return sum;
}

void Multitime::extract(std::vector<double>& t) const
{
  t.clear();

  for(int i=0; i<Multitime::tdimension; ++i) {
    if (!active[i]) continue;
    t.push_back(std::abs(v2[i] - v1[i]) + 0.5*(std::abs(v1[i]) + std::abs(v2[i])));
  }
}

namespace SYNARMOSMA {
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
    for(i=0; i<Multitime::tdimension; ++i) {
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
}

