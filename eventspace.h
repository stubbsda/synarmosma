#include "event.h"
#include "nexus.h"

#ifndef __eventspaceh
#define __eventspaceh

namespace SYNARMOSMA {
  typedef boost::adjacency_list<boost::vecS,boost::vecS,boost::bidirectionalS> Network;
  
  class Eventspace {
   private:
    Event* events;  
    int nevent;
    Network* causal_network;
    double L;

    void initiate_events();
    double perceptual_divergence(const double*,double,const double*,const double*) const; 
    double compute_distance(int,int) const;
    void perceptual_consistency();
   public:
    Eventspace();
    Eventspace(int);
    Eventspace(const Eventspace&);
    ~Eventspace();
    Eventspace& operator=(const Eventspace&);
    void compute_distances(std::vector<double>&) const;
    void write(const std::vector<double>&,double,const char*) const;
    void build_tessellation(Nexus*) const;
  };
}
#endif 
