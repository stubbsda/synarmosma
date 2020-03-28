#include "markov_chain.h"

using namespace SYNARMOSMA;

extern Random RND;

Markov_Chain::Markov_Chain(unsigned int n,double severity)
{
  if (n < 2) throw std::invalid_argument("The dimension of the Markov chain's state space must be greater than zero!");
  N = n;
  // The transition matrix will initially just be the identity matrix...
  transition_matrix = new Matrix<double>(N,true);
  mutate(severity);
}

Markov_Chain::~Markov_Chain()
{
  if (N > 0) delete transition_matrix;
}

void Markov_Chain::mutate(double severity)
{
  if (severity < -std::numeric_limits<double>::epsilon()) throw std::invalid_argument("The Markov_Chain::mutation argument must be non-negative!");
  // If the severity is zero there's nothing to do...
  if (severity < std::numeric_limits<double>::epsilon()) return;
  int r,n,m;
  bool done = false;
  double c,sigma,mtotal = 0.0;
  std::vector<double> vx;
  const double pfactor = severity/double(N);

  do {
    r = RND.irandom(N);
    transition_matrix->get_row(vx,r);
    n = RND.irandom(N);
    do {
      m = RND.irandom(N);
      if (m == n) continue;
      if (std::abs(vx[n] - vx[m]) > std::numeric_limits<double>::epsilon()) break;
    } while(true);
    sigma = pfactor*RND.drandom();
    if (vx[n] > vx[m]) {
      c = sigma*vx[n];
      c = std::min(c,1.0 - vx[m]); 
      vx[n] -= c;
      vx[m] += c;
    }
    else {
      c = sigma*vx[m];
      c = std::min(c,1.0 - vx[n]); 
      vx[m] -= c;
      vx[n] += c;
    }
    mtotal += c;
    if (mtotal > pfactor) done = true;
  } while(!done);

  assert(consistent());
}

bool Markov_Chain::consistent() const
{
  if (N < 2) return false;
  if (transition_matrix->get_nrow() != N || transition_matrix->get_ncolumn() != N) return false;
  unsigned int i,j;
  double sum;
  std::vector<double> vx;

  for(i=0; i<N; ++i) {
    transition_matrix->get_row(vx,i);
    sum = 0.0;
    for(j=0; j<N; ++j) {
      if (vx[j] < -std::numeric_limits<double>::epsilon()) return false;
      sum += vx[j];
    }
    if (std::abs(sum - 1.0) > std::numeric_limits<double>::epsilon()) return false;
  }
  return true;
}

void Markov_Chain::multiply(const std::vector<double>& cstate,std::vector<double>& output) const
{
  unsigned int i,j;
  double sum;
  std::vector<double> vx;

  output.clear();
  for(i=0; i<N; ++i) {
    transition_matrix->get_column(vx,i);
    sum = 0.0;
    for(j=0; j<N; ++j) {
      sum += cstate[j]*vx[j];
    }
    output.push_back(sum);
  }
}

void Markov_Chain::set_transition(const std::vector<double>& vx,unsigned int n)
{
  if (n >= N) throw std::invalid_argument("Illegal state value in Markov_Chain::set_transition!");
  if (vx.size() != N) throw std::invalid_argument("Illegal row vector length in Markov_Chain::set_transition!");
  transition_matrix->set_row(vx,n);
  assert(consistent());
}

void Markov_Chain::get_state(const std::vector<double>& cstate,std::vector<double>& output,int n) const
{
  if (cstate.size() != N) throw std::invalid_argument("Illegal state vector length in Markov_Chain::get_state!");
  if (n < 1) throw std::invalid_argument("The number of steps in Markov_Chain::get_state must be greater than zero!");  
  if (n == 1) {
    multiply(cstate,output);
  }
  else {
    unsigned int i,j;
    double sum;
    std::vector<double> vx;
    Matrix<double> tpower = transition_matrix->pow(n);

    output.clear();
    for(i=0; i<N; ++i) {
      tpower.get_column(vx,i);
      sum = 0.0;
      for(j=0; j<N; ++j) {
        sum += cstate[j]*vx[j];
      }
      output.push_back(sum);
    }    
  }
}

int Markov_Chain::get_state(const std::vector<double>& cstate,int n) const
{
  unsigned int i;
  int output = -1;
  double sum = 0.0,alpha = RND.drandom();
  std::vector<double> pvector;

  get_state(cstate,pvector,n);

  for(i=0; i<N; ++i) {
    sum += pvector[i];
    if (alpha < sum) {
      output = i;
      break;
    }
  }
  if (output == -1) output = N - 1;

  return output;
}

int Markov_Chain::get_state(unsigned int n) const
{
  if (n >= N) throw std::invalid_argument("Illegal state value in Markov_Chain::get_state!");
  unsigned int i;
  int output = -1;
  double sum = 0.0,alpha = RND.drandom();
  std::vector<double> pvector;

  transition_matrix->get_row(pvector,n);

  for(i=0; i<N; ++i) {
    sum += pvector[i];
    if (alpha < sum) {
      output = i;
      break;
    }
  }
  if (output == -1) output = N - 1;

  return output;
}