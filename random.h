#include "global.h"

#ifndef _randomh
#define _randomh

namespace SYNARMOSMA {
  class Random {
   private:
    unsigned int s;
    base_generator_type BGT;

    boost::math::beta_distribution<>* root_beta;

    boost::bernoulli_distribution<>* brn;
    boost::variate_generator<base_generator_type&,boost::bernoulli_distribution<> >* vbrn;

    boost::poisson_distribution<>* poisson;
    boost::variate_generator<base_generator_type&,boost::poisson_distribution<> >* vpoisson;

    boost::uniform_real<>* uniform;
    boost::variate_generator<base_generator_type&,boost::uniform_real<> >* VRG;

    boost::normal_distribution<>* gaussian;
    boost::variate_generator<base_generator_type&,boost::normal_distribution<> >* NRG;

    bool brn_allocated;
    bool beta_allocated;
    bool poisson_allocated;
   public:
    Random();
    ~Random();
    inline void set_seed(unsigned int x) {s = x; BGT.seed(s);};
    inline void increment_seed() {s++; BGT.seed(s);};
    inline void decrement_seed() {s--; BGT.seed(s);};
    inline unsigned int get_seed() const {return s;};
    void initialize_beta(double,double);
    void initialize_bernoulli(double);
    void initialize_poisson(double);
    double drandom();
    double drandom(double,double);
    double beta_variate();
    bool bernoulli_variate();
    bool poisson_variate();
    double nrandom();
    double nrandom(double);
    double nrandom(double,double);
    int irandom(int);
    int irandom(int,int);
    int irandom(const std::set<int>&);
    unsigned int irandom(const std::set<unsigned int>&);
    int irandom(const std::set<int>&,const std::set<int>&);
    int irandom(const std::vector<int>&);
    int irandom(int,const std::vector<int>&);
    void shuffle(std::vector<int>&,int);
    void generate_random_vector(std::vector<double>&,int,double,double);
    void generate_random_vector(std::vector<std::complex<double> >&,int,double,double,bool = false);
  };
}
#endif