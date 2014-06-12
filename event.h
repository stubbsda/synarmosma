#include "proposition.h"
#include "multitime.h"
#include "geometry.h"

#ifndef _eventh
#define _eventh

class Event {
 private:
  int timestep;
  int colour;
  int local_dimension;
  double energy;
  Proposition observation;
  std::set<int> future,past;  
  std::set<int> neighbours;
  std::vector<std::string> entourage;
  std::vector<double> space;  
  Multitime proper_time;

  double norm() const;
  void allocate();
  void initialize();
 public:
  Event();
  Event(const Event&);
  Event(const char*);
  ~Event();  
  Event& operator=(const Event&);
  friend class Eventspace;
  friend double compute_distance(const Event&,const Event&);
};

// A class which implements a Leibnizian concept of spacetime 
/*
class Event_Graph: public Graph {
 private:    
  std::vector<unsigned int> events;     

  void perceptual_consistency();
  kind perceptual_divergence(double*,double,double*,double**) const; 
  void allocate();
  void initialize();
  
 public:
  Event_Graph();
  Event_Graph(const char*);
  Event_Graph(const Event_Graph&);
  Event_Graph & operator=(const Event_Graph&);
  ~Event_Graph();
  void initiate_events();
};
*/
#endif 

