#include "random.h"

using namespace SYNARMOSMA;

Random::Random()
{
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
  if (brn_allocated) {
    delete brn;
    delete vbrn;
  }
  else {
    brn_allocated = true;
  }
  brn = new boost::bernoulli_distribution<>(p);
  vbrn = new boost::variate_generator<base_generator_type&,boost::bernoulli_distribution<> >(BGT,*brn);
}

void Random::initialize_poisson(double lambda)
{
  if (poisson_allocated) {
    delete poisson;
    delete vpoisson;
  }
  else {
    poisson_allocated = true;
  }
  poisson = new boost::poisson_distribution<>(lambda);
  vpoisson = new boost::variate_generator<base_generator_type&,boost::poisson_distribution<> >(BGT,*poisson);
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
  if (y < x) throw std::invalid_argument("The arguments to Random::drandom(a,b) must satisfy a < b!");
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
  if (sigma < std::numeric_limits<double>::epsilon()) throw std::invalid_argument("The variance of the Gaussian distribution must be greater than zero!");
  double z = sigma*nrandom() + mu;
  return z;
}

double Random::nrandom(double sigma)
{
  if (sigma < std::numeric_limits<double>::epsilon()) throw std::invalid_argument("The variance of the Gaussian distribution must be greater than zero!");
  double z = nrandom(0.0,sigma);
  return z;
}

int Random::irandom(int n)
{
  if (n < 1) throw std::invalid_argument("The argument of Random::irandom must be greater than zero!");
  int output = int(n*drandom());
  return output;
}

int Random::irandom(int x,int y)
{
  if (x == y) return x;
  int output = x + irandom(y-x);
  return output;
}

unsigned int Random::irandom(const std::set<unsigned int>& S)
{
  if (S.empty()) throw std::invalid_argument("The argument of Random::irandom must not be an empty set!");
  unsigned int output = 0;
  int n = irandom(S.size());
  int k = 0;
  std::set<unsigned int>::const_iterator it;

  for(it=S.begin(); it!=S.end(); ++it) {
    if (k == n) {
      output = *it;
      break;
    }
    k++;
  }
  return output;
}

int Random::irandom(const std::set<int>& S)
{
  if (S.empty()) throw std::invalid_argument("The argument of Random::irandom must not be an empty set!");
  int output = 0;
  unsigned int k = 0,n = irandom(S.size());
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
  if (S.empty()) throw std::invalid_argument("The argument of Random::irandom must not be an empty set!");
  int output = 0;

  do {
    output = irandom(S);
    if (used.count(output) == 0) break;
  } while(true);
  return output;
}

int Random::irandom(const std::vector<int>& used)
{
  if (std::count(used.begin(),used.end(),0) == 0) throw std::invalid_argument("There must be at least one zero element in the argument of Random:irandom!"); 
  unsigned int l = used.size();
  int output = 0;

  do {
    output = irandom(l);
    if (used[output] == 0) break;
  } while(true);
  return output;
}

int Random::irandom(int ulimit,const std::vector<int>& used)
{
  int output = 0;

  do {
    output = irandom(ulimit);
    if (std::count(used.begin(),used.end(),output) == 0) break;
  } while(true);
  return output;
}

void Random::shuffle(std::vector<int>& index,int n)
{
  if (n < 1) throw std::invalid_argument("The number of elements to shuffle must be greater than zero!");
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

void Random::generate_random_vector(std::vector<double>& x,int n,double a,double b) 
{
  if (n < 1) throw std::invalid_argument("The length of the random vector must be greater than zero!");

  x.clear();

  for(int i=0; i<n; ++i) {
    x.push_back(drandom(a,b));
  }
}

void Random::generate_random_vector(std::vector<std::complex<double> >& x,int n,double a,double b,bool nz_imaginary) 
{
  if (n < 1) throw std::invalid_argument("The length of the random vector must be greater than zero!");

  std::complex<double> alpha; 

  x.clear();

  if (nz_imaginary) {
    for(int i=0; i<n; ++i) {
      alpha = std::complex<double>(drandom(a,b),drandom(a,b));
      x.push_back(alpha);
    }
  }
  else {
    for(int i=0; i<n; ++i) {
      alpha = std::complex<double>(drandom(a,b),0.0);
      x.push_back(alpha);
    }
  }
}

