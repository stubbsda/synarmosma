#include "eventspace.h"

Eventspace::Eventspace()
{
  nevent = 20;
  events = new Event[nevent];
}

Eventspace::Eventspace(int n)
{
  nevent = n;
  events = new Event[nevent];
}

Eventspace::Eventspace(const Eventspace& source)
{
  nevent = source.nevent;
  L = source.L;
  // Is this a deep copy?
  causal_network = source.causal_network;
  events = new Event[nevent];
  for(int i=0; i<nevent; ++i) {
    events[i] = source.events[i];
  }
}

Eventspace& Eventspace::operator=(const Eventspace& source)
{
  if (this != &source) {
    nevent = source.nevent;
    L = source.L;
    // Is this a deep copy?
    causal_network = source.causal_network;
    events = new Event[nevent];
    for(int i=0; i<nevent; ++i) {
      events[i] = source.events[i];
    }
  }
  return *this;
}

Eventspace::~Eventspace()
{
  delete[] events;
}

void Eventspace::build_tessellation(Nexus* output) const
{
  int i,j;
  std::vector<double> coordinates;

  for(i=0; i<nevent; ++i) {
    events[i].proper_time.extract(coordinates);
    for(j=0; j<events[i].local_dimension; ++j) {
      coordinates.push_back(events[i].space[j]);
    }
  }
  output->clear();
  // What now to build the Nexus instance?

}

void Eventspace::compute_distances(std::vector<double>& distances) const
{
  int i,j,k,l;
  double delta;
  std::vector<double> t1,t2;

  for(i=0; i<nevent; ++i) {
    events[i].proper_time.extract(t1);
    for(j=1+i; j<nevent; ++j) {
      events[j].proper_time.extract(t2);
      delta = 0.0;
      for(k=0; k<tdimension; ++k) {
        delta += (t1[k]-t2[k])*(t1[k]-t2[k]);
      }
      l = (events[i].local_dimension < events[j].local_dimension) ? events[i].local_dimension : events[j].local_dimension;
      for(k=0; k<l; ++k) {
        delta += (events[i].space[k] - events[j].space[k])*(events[i].space[k] - events[j].space[k]);
      }
      delta = std::sqrt(delta);
      distances.push_back(delta);
    }
  }
}

void Eventspace::write(const std::vector<double>& distances,double delta,const char* filename) const
{
  std::ofstream s(filename,std::ios_base::out | std::ios_base::trunc | std::ios_base::binary);
  int i,j,in1 = 0,neg1 = -1;
  float x;

  s.write((char*)(&nevent),sizeof(int));
  for(i=0; i<nevent; ++i) {
    // What should be done if the the number of spatial coordinates isn't equal to three?
    for(j=0; j<events[i].local_dimension; ++j) {
      x = float(events[i].space[j]);
      s.write((char*)(&x),sizeof(float));
    }
  }

  for(i=0; i<nevent; ++i) {
    for(j=i+1; j<nevent; ++j) {
      if (distances[in1] < delta) s.write((char*)(&j),sizeof(int));
      in1 += 1;
    }
    s.write((char*)(&neg1),sizeof(int));

  }
  s.close();
}

void Eventspace::perceptual_consistency()
{

}

double Eventspace::perceptual_divergence(const double* raxis,double theta,const double* translation,const double* observed) const
{
  // The method assumes that
  // a) raxis is a unit vector 
  // b) 0 <= theta < 2*pi
  // c) dimension = 3
  int i,j;
  double ct,st,d1,d2,xt,q0,temp[3],delta,sum = 0.0;
  const double l2 = 2.0*L;
  const double zero = 0.0;

  ct = std::cos(theta);
  st = std::sin(theta);
#ifndef USE_MPI
#ifdef PARALLEL
#pragma omp parallel for default(shared) private(i,j,xt,temp,delta,d1,d2,q0) reduction(+:sum)
#endif
#endif
  for(i=0; i<nevent; ++i) {
    // The angular motion:
    // We calculate the cross product...
    temp[0] = events[i].space[1]*raxis[2]-events[i].space[2]*raxis[1];
    temp[1] = events[i].space[2]*raxis[0]-events[i].space[0]*raxis[0];
    temp[2] = events[i].space[0]*raxis[1]-events[i].space[1]*raxis[0];
    // And now the scalar product
    q0 = events[i].space[0]*raxis[0] + events[i].space[1]*raxis[1] + events[i].space[2]*raxis[2];
    delta = zero;
    for(j=0; j<3; ++j) {
      xt = q0*raxis[j] + st*temp[j] + ct*(events[i].space[j]-q0*raxis[j]) + translation[j];
      d1 = xt - observed[3*i+j];
      if (d1 < zero) d1 = -d1;
      d2 = l2 - d1;
      if (d1 > d2) d1 = d2;
      delta += d1*d1;
    }
    sum += std::sqrt(delta);
  }
  return (sum/double(nevent));
}

void Eventspace::initiate_events()
{

}
