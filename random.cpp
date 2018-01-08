#include "global.h"

using namespace SYNARMOSMA;

Random::Random()
{
  brn_allocated = false;
  beta_allocated = false;
  poisson_allocated = false;
  uniform = new boost::uniform_real<>(0,1);
  gaussian = new boost::normal_distribution<>;
  VRG = new boost::variate_generator<base_generator_type&,boost::uniform_real<> >(BGT,*uniform);
  NRG = new boost::variate_generator<base_generator_type&,boost::normal_distribution<> >(BGT,*gaussian);
}

Random::~Random()
{
  delete uniform;
  delete gaussian;
  delete VRG;
  delete NRG;
  if (poisson_allocated) {
    delete poisson;
    delete vpoisson;
  }
  if (brn_allocated) {
    delete brn;
    delete vbrn;
  }
  if (beta_allocated) delete root_beta;
}

void Random::initialize_bernoulli(double p)
{
  brn = new boost::bernoulli_distribution<>(p);
  vbrn = new boost::variate_generator<base_generator_type&,boost::bernoulli_distribution<> >(BGT,*brn);
  brn_allocated = true;
}

void Random::initialize_poisson(double lambda)
{
  poisson = new boost::poisson_distribution<>(lambda);
  vpoisson = new boost::variate_generator<base_generator_type&,boost::poisson_distribution<> >(BGT,*poisson);
  poisson_allocated = true;
}

void Random::initialize_beta(double a,double b)
{
  if (beta_allocated) {
    delete root_beta;
  }
  else {
    beta_allocated = true;
  }
  root_beta = new boost::math::beta_distribution<>(a,b);
}

double Random::drandom()
{
  // Returns a uniform random variate on [0,1)
  return (*VRG)();
}

double Random::drandom(double x,double y)
{
  double alpha = x + drandom()*(y - x);
  return alpha;
}

double Random::beta_variate()
{
  return boost::math::quantile(*root_beta,drandom());
}

bool Random::bernoulli_variate()
{
  return (*vbrn)();
}

bool Random::poisson_variate()
{
  return (*vpoisson)();
}

double Random::nrandom()
{
  return (*NRG)();
}

double Random::nrandom(double mu,double sigma)
{
  double z = sigma*nrandom() + mu;
  return z;
}

double Random::nrandom(double sigma)
{
  double z = nrandom(0.0,sigma);
  return z;
}

int Random::irandom(int n)
{
  int output = int(n*drandom());
  return output;
}

int Random::irandom(int x,int y)
{
  int output = irandom(y-x);
  return (x + output);
}

int Random::irandom(const std::set<int>& S)
{
  int output = 0;
  int n = irandom(S.size());
  int k = 0;
  std::set<int>::const_iterator it;

  for(it=S.begin(); it!=S.end(); ++it) {
    if (k == n) {
      output = *it;
      break;
    }
    k++;
  }
  return output;
}

int Random::irandom(const std::set<int>& S,const std::set<int>& used)
{
  int output = 0;
  std::set<int>::const_iterator it;

  do {
    output = irandom(S);
    it = used.find(output);
    if (it == used.end()) break;
  } while(true);
  return output;
}

int Random::irandom(const std::vector<int>& used)
{
  int i,output,l = (signed) used.size();
  bool good = false;

  // But what if used[i] = 1 for all i?
  for(i=0; i<l; ++i) {
    if (used[i] == 0) {
      good = true;
      break;
    }
  }
  if (!good) return irandom(l);
  do {
    output = irandom(l);
    if (used[output] == 0) break;
  } while(true);
  return output;
}

int Random::irandom(int ulimit,const std::vector<int>& used)
{
  int i,output = 0,l = (signed) used.size();
  bool found;
  do {
    output = irandom(ulimit);
    found = false;
    for(i=0; i<l; ++i) {
      if (used[i] == output) {
        found = true;
        break;
      }
    }
    if (!found) break;
  } while(true);
  return output;
}

void Random::shuffle(std::vector<int>& index,int n)
{
  int i,r,t;
  
  index.clear();
  for(i=0; i<n; ++i) {
    index.push_back(i);
  }

  for(i=0; i<n-1; ++i) {
    r = i + irandom(0,n-i);
    t = index[i];
    index[i] = index[r];
    index[r] = t;
  }
}
