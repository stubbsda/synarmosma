#include "homotopy.h"

Homotopy::Homotopy()
{
  unsigned int i,n = 15;
  for(i=0; i<n; ++i) {
    Group g;
    sequence.push_back(g);
  }
  compute_fitness();
}

Homotopy::Homotopy(const Homotopy& source)
{
  sequence = source.sequence;
  fitness = source.fitness;
}

Homotopy::Homotopy(unsigned int n)
{
  unsigned int i;
  for(i=0; i<n; ++i) {
    Group g;
    sequence.push_back(g);
  }
  compute_fitness();
}

Homotopy::~Homotopy()
{

}

void Homotopy::compute_fitness()
{
  unsigned int i;
  float temp,sum = 0.0;
  for(i=0; i<sequence.size(); ++i) {
    temp = std::exp(-std::pow(float(sequence[i].ngenerator)-3.0,2))+ std::pow(std::sin(float(sequence[i].relations.size())),2);
    sum += temp;
  }
  fitness = sum;
}

float Homotopy::get_fitness() const
{
  return fitness;
}

void Homotopy::mutate()
{
  unsigned int n = irandom(sequence.size());
  Group g;
  sequence[n] = g;
  compute_fitness();
}

void Homotopy::write()
{
  std::cout << "This homotopy sequence has " << sequence.size() << " groups and a fitness of " << fitness << std::endl;
}

Homotopy& Homotopy::operator =(const Homotopy& source)
{
  if (this == &source) return *this;

  sequence = source.sequence;
  fitness = source.fitness;

  return *this;
}

Homotopy operator +(const Homotopy& h1,const Homotopy& h2)
{
  if (h1.sequence.size() != h2.sequence.size()) {
    std::cerr << "These two homotopy sequences cannot be added: they don't have the same length!" << std::endl;
    std::exit(1);
  }
  unsigned int i,bisection;

  Homotopy output;
  bisection = 1 + irandom(h1.sequence.size()-1);
  output.sequence.clear();
  for(i=0; i<bisection; ++i) {
    output.sequence.push_back(h1.sequence[i]);
  }
  for(i=bisection; i<h1.sequence.size(); ++i) {
    output.sequence.push_back(h2.sequence[i]);
  }

  return output;
}

