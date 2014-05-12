#include "event.h"

extern Random RND;

Event::Event()
{
  allocate();
  initialize();
}

Event::Event(const char* filename)
{
  allocate();
  initialize();
}

Event::Event(const Event& source)
{
#ifndef LEIBNIZ
  proper_time = source.proper_time;
  space = source.space;
#endif
  energy = source.energy;
  neighbours = source.neighbours;
  entourage = source.entourage;
  observation = source.observation;
  timestep = source.timestep;
  local_dimension = source.local_dimension;
  colour = source.colour;
  past = source.past;
  future = source.future;
}

Event& Event::operator=(const Event& source)
{
  if (this != &source) {
#ifndef LEIBNIZ
    proper_time = source.proper_time;
    space = source.space;
#endif
    energy = source.energy;
    neighbours = source.neighbours;
    entourage = source.entourage;
    observation = source.observation;
    timestep = source.timestep;
    local_dimension = source.local_dimension;
    colour = source.colour;
    past = source.past;
    future = source.future;
  }
  return *this;
}

Event::~Event()
{

}

void Event::allocate()
{


}

void Event::initialize()
{
  local_dimension = 0;
  timestep = 0;
  colour = 2;
#ifndef LEIBNIZ
  proper_time.initialize(RND.nrandom(1.5));
  for(int i=0; i<Geometry::background_dimension; ++i) {
    space.push_back(-5.0 + RND.nrandom(10.0));
  }
#endif
}

double Event::norm() const
{
  double output = 0.0;
#ifndef LEIBNIZ
  double rl = 0.0;
  for(int i=0; i<Geometry::background_dimension; ++i) {
    rl += std::abs(space[i])*std::abs(space[i]);
  }
  rl = std::sqrt(rl);
  output = proper_time.norm() + rl;
#endif
  return output;
}

double compute_distance(const Event& e1,const Event& e2)
{
  double output = 0.0;
#ifdef LEIBNIZ

#else
  int i,n = (e1.space.size() <= e2.space.size()) ? e1.space.size() : e2.space.size();
  double t1,t2,rl = 0.0;
  e1.proper_time.extract(&t1);
  e2.proper_time.extract(&t2);
  for(i=0; i<n; ++i) {
    rl += std::abs(e1.space[i] - e2.space[i])*std::abs(e1.space[i] - e2.space[i]);
  }
  output = (t1 - t2);
  output = output*output;
  output += rl;
  output = std::sqrt(output);
#endif
  return output;
}








