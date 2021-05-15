#include "functional_equation.h"

using namespace SYNARMOSMA;

Functional_Equation::Functional_Equation()
{

}

Functional_Equation::Functional_Equation(const std::string& filename)
{
  unsigned int i;
  std::string line,store;
  std::vector<std::string> alpha,beta,exponents;
  std::vector<unsigned int> bk;

  std::ifstream s;
  s.exceptions(std::ifstream::badbit);
  try {
    s.open(filename,std::ios::in);

    // Loop through all lines in the parameter file
    while(std::getline(s,line)) {
      // If it's an empty line, continue
      if (line.empty()) continue;
      // If the line begins with a #, ignore it
      if (line[0] == '#') continue;
      // Find the position of the colons
      for(i=0; i<line.length(); ++i) {
        if (line[i] == ':') bk.push_back(i);
      }
      if (bk.size() != 2) throw std::runtime_error("Error parsing the functional equation input file!");
      store = line.substr(0,bk[0]);
      trim(store);
      alpha.push_back(store);
      store = line.substr(bk[0]+1,bk[1]-1);
      trim(store);
      beta.push_back(store);
      store = line.substr(bk[1]+1,line.length());
      trim(store);
      exponents.push_back(store);
      bk.clear();
    }
  }
  catch (const std::ifstream::failure& e) {
    std::cout << "Error in opening or reading the " << filename << " file!" << std::endl;
  }
  s.close();

  parse_equation(alpha,beta,exponents);
}

void Functional_Equation::parse_polynomial(const std::string& polynomial,std::vector<Rational>& coefficients) const
{
  // Note that this method assumes that a polynomial like 3x^4 - 11/2 x^2 + 7x - 2 is written as 
  // (-2,7,-11/2,0,3), i.e. in ascending order of degree.
  unsigned int i;
  int p,q;
  char c;
  bool top = true;
  std::string numerator,denominator;

  coefficients.clear();

  for(i=0; i<polynomial.length(); ++i) {
    c = polynomial[i];
    if (c == '(') {
      top = true;
      continue;
    }
    if (c == ',' || c == ')') {
      if (denominator.empty()) denominator = "1";
      p = boost::lexical_cast<signed int>(numerator);
      q = boost::lexical_cast<signed int>(denominator);
      coefficients.push_back(Rational(p,q));
      denominator.clear();
      numerator.clear();
      top = true;
      continue;
    }
    if (c == '/') {
      top = false;
      continue;
    }
    if (top) {
      numerator.push_back(c);
    }
    else {
      denominator.push_back(c);
    }
  }
}

void Functional_Equation::parse_equation(const std::vector<std::string>& alpha,const std::vector<std::string>& beta,const std::vector<std::string>& exponents)
{  
  unsigned int i,degree;
  Polynomial<Rational> p1,p2;
  std::vector<Rational> vx;

  for(i=0; i<alpha.size(); ++i) {
    degree = boost::lexical_cast<unsigned int>(exponents[i]);
    parse_polynomial(alpha[i],vx);
    p1 = Polynomial<Rational>(vx);
    if (degree == 0) {
      constant = p1;
      continue;
    }
    parse_polynomial(beta[i],vx);
    p2 = Polynomial<Rational>(vx);
    std::tuple<Polynomial<Rational>,Polynomial<Rational>,unsigned int> trio(p1,p2,degree);
    terms.push_back(trio);
  }
  simplify();
  if (!consistent()) throw std::invalid_argument("The functional equation input is not consistent with this class!");
}

Functional_Equation::Functional_Equation(const Functional_Equation& source)
{
  linear = source.linear;
  homogeneous = source.homogeneous;
  terms = source.terms;
  constant = source.constant;
}

Functional_Equation& Functional_Equation::operator =(const Functional_Equation& source)
{
  if (this == &source) return *this;

  linear = source.linear;
  homogeneous = source.homogeneous;
  terms = source.terms;
  constant = source.constant;

  return *this;
}

Functional_Equation::~Functional_Equation()
{

}

void Functional_Equation::clear()
{
  terms.clear();
  constant.clear();
  linear = false;
  homogeneous = false;
}

int Functional_Equation::serialize(std::ofstream& s) const
{
  unsigned int i,j,n;
  int count = 0;
  Polynomial<Rational> p;
  
  s.write((char*)(&linear),sizeof(bool)); count += sizeof(bool);
  s.write((char*)(&homogeneous),sizeof(bool)); count += sizeof(bool);
  n = terms.size();
  s.write((char*)(&n),sizeof(int)); count += sizeof(int);
  for(i=0; i<n; ++i) {
    j = std::get<2>(terms[i]);
    s.write((char*)(&j),sizeof(int)); count += sizeof(int);
    p = std::get<0>(terms[i]);
    count += p.serialize(s);
    p = std::get<1>(terms[i]);
    count += p.serialize(s);
  }
  count += constant.serialize(s);
  
  return count;
}

int Functional_Equation::deserialize(std::ifstream& s)
{
  unsigned int i,j,n;
  int count = 0;
  Polynomial<Rational> p1,p2;

  clear();

  s.read((char*)(&linear),sizeof(bool)); count += sizeof(bool);
  s.read((char*)(&homogeneous),sizeof(bool)); count += sizeof(bool);
  s.read((char*)(&n),sizeof(int)); count += sizeof(int);
  for(i=0; i<n; ++i) {
    s.read((char*)(&j),sizeof(int)); count += sizeof(int);
    count += p1.deserialize(s);
    count += p2.deserialize(s);
    terms.push_back(std::tuple<Polynomial<Rational>,Polynomial<Rational>,unsigned int>(p1,p2,j));
  }
  count += constant.deserialize(s);

  return count;
}

bool Functional_Equation::consistent() const
{
  unsigned int i,d,n = terms.size();
  std::set<unsigned int> dset;

  // Begin by checking that all of the terms are distinct in terms of their degree:
  for(i=0; i<n; ++i) {
    d = std::get<2>(terms[i]);
    if (dset.count(d) > 0) return false;
    dset.insert(d);
  }
  return true;
}

bool Functional_Equation::simplify()
{
  unsigned int i,d = 0,n = terms.size();
  bool output = false;
  Polynomial<Rational> q;
  std::vector<std::tuple<Polynomial<Rational>,Polynomial<Rational>,unsigned int> > nterms;

  if (!homogeneous) {
    if (constant.is_null()) {
      homogeneous = true;
      output = true;
    }
  }
  for(i=0; i<n; ++i) {
    q = std::get<0>(terms[i]);
    if (q.is_null()) continue;
    if (std::get<2>(terms[i]) > d) d = std::get<2>(terms[i]);
    nterms.push_back(terms[i]);
  }
  if (d == 1) {
    if (!linear) {
      linear = true;
      output = true;
    }
  }
  if (nterms.size() == terms.size()) return output;
  terms = nterms;
  return true;  
}

bool Functional_Equation::reduce(unsigned int p,Variety<unsigned int>* output) const
{
  if (!NTL::ProbPrime(p)) throw std::invalid_argument("Functional equation must be reduced over a prime!");

  unsigned int i,j,degree = terms.size();
  bool failed;
  NTL::ZZ z;
  std::pair<unsigned int,unsigned int> duo;
  std::vector<unsigned int> solutions;
  Rational q,a,gamma;
  Monomial<unsigned int> t;
  std::vector<Rational> cvector,dvector;

  output->initialize(p,p);

  for(i=0; i<p; ++i) {
    q = Rational(i,1);
    for(j=0; j<degree; ++j) {
      a = std::get<0>(terms[j]).evaluate(q);
      cvector.push_back(a);
      a = std::get<1>(terms[j]).evaluate(q);
      dvector.push_back(a);
    }
    gamma = constant.evaluate(q);
    // If any of the delay vector elements has a denominator 
    // which is divisible by the field characteristic, this 
    // method must return with the information that no such 
    // reduction is possible.
    failed = false;
    for(j=0; j<degree; ++j) {
      z = dvector[j].get_denominator();
      if (z % p == 0) {
        failed = true;
        break;
      }
    }
    if (failed) return false;
    z = gamma.get_denominator();
    for(j=0; j<degree; ++j) {
      z *= cvector[j].get_denominator();
    }
    a = Rational(z,NTL::ZZ(1));
    for(j=0; j<degree; ++j) {
      q = a*cvector[j];
      t.coefficient = q.reduce(p);
      duo.first = dvector[j].reduce(p);
      duo.second = std::get<2>(terms[j]);
      t.exponents.push_back(duo);
      output->add_term(i,t);
    }
    q = a*gamma;
    output->set_constant(i,q.reduce(p)); 
  }

  output->compute_properties();
  return true;
}
