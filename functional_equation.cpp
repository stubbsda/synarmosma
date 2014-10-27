#include "functional_equation.h"

using namespace SYNARMOSMA;

extern Random RND;

template<class kind>
Functional_Equation<kind>::Functional_Equation()
{
  unsigned int n = 5 + RND.irandom(6);

  linear = false;
  homogeneous = false;
  float alpha = float(RND.drandom());
  if (alpha < 0.5) linear = true;
  alpha = float(RND.drandom());
  if (alpha < 0.5) homogeneous = true;
  initialize(n);
}

template<class kind>
Functional_Equation<kind>::Functional_Equation(unsigned int n)
{
  initialize(n);
}

template<class kind>
Functional_Equation<kind>::Functional_Equation(const char* filename)
{
  unsigned int i;
  bool first = true;
  std::string line,store;
  std::vector<std::string> alpha,beta,exponents;
  std::vector<unsigned int> bk;

  std::ifstream s(filename,std::ios_base::in);
  if (!s.is_open()) {
    // File doesn't exist, print an error message and die
    std::cerr << "The file " << filename << " cannot be found!" << std::endl;
    std::exit(1);
  }
  // Loop through all lines in the parameter file
  while(std::getline(s,line)) {
    // If it's an empty line, continue
    if (line.empty()) continue;
    // If the line begins with a #, ignore it
    if (line[0] == '#') continue;
    if (first) {
      // We need to determine the field of this equation 
      trim(line);
      //if (linie[0] == 'Q') koerper = 1;
      first = false;
      continue;
    }
    // Find the position of the colons
    for(i=0; i<line.length(); ++i) {
      if (line[i] == ':') bk.push_back(i);
    }
    if (bk.size() != 2) {
      std::cerr << "The file format is incorrect at line " << line << std::endl;
      std::exit(1);
    }
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
  s.close();

  analyze_file(alpha,beta,exponents);
}

namespace SYNARMOSMA {
template<>
void Functional_Equation<Rational>::analyze_file(std::vector<std::string>& alpha,std::vector<std::string>& beta,std::vector<std::string>& exponents)
{
  Polynomial<Rational> p1,p2;
  std::vector<Rational> vx;
  Rational coefficient;
  unsigned int i,j,degree;
  char c;
  std::string s1,s2;
  bool numerator;

  for(i=0; i<alpha.size(); ++i) {
    degree = boost::lexical_cast<unsigned int>(exponents[i]);
    numerator = true;
    for(j=0; j<alpha[i].length(); ++j) {
      c = alpha[i][j];
      if (c == '(') continue;
      if (c == ',' || c == ')') {
        coefficient = Rational(boost::lexical_cast<signed int>(s1),boost::lexical_cast<signed int>(s2));
        vx.push_back(coefficient);
        s1.clear();
        s2.clear();
        continue;
      }
      if (c == '/') {
        numerator = false;
        continue;
      }
      if (numerator) {
        s1.push_back(c);
      }
      else {
        s2.push_back(c);
      }
    }
    p1 = Polynomial<Rational>(vx);

    vx.clear();
    s1.clear();
    s2.clear();
    if (degree == 0) {
      remainder = p1;
      continue;
    }
    numerator = true;
    for(j=0; j<beta[i].length(); ++j) {
      c = beta[i][j];
      if (c == '(') continue;
      if (c == ',' || c == ')') {
        coefficient = Rational(boost::lexical_cast<signed int>(s1),boost::lexical_cast<signed int>(s2));
        vx.push_back(coefficient);
        s1.clear();
        s2.clear();
        continue;
      }
      if (c == '/') {
        numerator = false;
        continue;
      }
      if (numerator) {
        s1.push_back(c);
      }
      else {
        s2.push_back(c);
      }
    }
    p2 = Polynomial<Rational>(vx);

    boost::tuple<Polynomial<Rational>,Polynomial<Rational>,unsigned int> trio(p1,p2,degree);
    terms.push_back(trio);
    vx.clear();
    s1.clear();
    s2.clear();
  }
}
}

template<class kind>
void Functional_Equation<kind>::analyze_file(std::vector<std::string>& alpha,std::vector<std::string>& beta,std::vector<std::string>& exponents)
{
  unsigned int i,j,degree;
  char c;
  std::string store;
  Polynomial<signed int> p1,p2;
  std::vector<signed int> vx;

  for(i=0; i<alpha.size(); ++i) {
    degree = boost::lexical_cast<unsigned int>(exponents[i]);
    for(j=0; j<alpha[i].length(); ++j) {
      c = alpha[i][j];
      if (c == '(') continue;
      if (c == ',' || c == ')') {
        vx.push_back(boost::lexical_cast<signed int>(store));
        store.clear();
        continue;
      }
      store.push_back(c);
    }
    p1 = Polynomial<signed int>(vx);

    vx.clear();
    store.clear();
    if (degree == 0) {
      remainder = p1;
      continue;
    }
    for(j=0; j<beta[i].length(); ++j) {
      c = beta[i][j];
      if (c == '(') continue;
      if (c == ',' || c == ')') {
        vx.push_back(boost::lexical_cast<signed int>(store));
        store.clear();
        continue;
      }
      store.push_back(c);
    }
    p2 = Polynomial<signed int>(vx);

    vx.clear();
    store.clear();
    boost::tuple<Polynomial<signed int>,Polynomial<signed int>,unsigned int> trio(p1,p2,degree);
    terms.push_back(trio);
  }
}

template<class kind>
Functional_Equation<kind>::Functional_Equation(const Functional_Equation& source)
{
  terms = source.terms;
}

template<class kind>
Functional_Equation<kind>& Functional_Equation<kind>::operator =(const Functional_Equation<kind>& source)
{
  if (this == &source) return *this;
  unsigned int i;
  linear = source.linear;
  homogeneous = source.homogeneous;
  for(i=0; i<source.terms.size(); ++i) {
    terms.push_back(source.terms[i]);
  }
  return *this;
}

template<class kind>
Functional_Equation<kind>::~Functional_Equation()
{

}

template<class kind>
void Functional_Equation<kind>::initialize(unsigned int n)
{
  unsigned int i,m = 1;
  for(i=1; i<n; ++i) {
    Polynomial<kind> alpha,beta;
    if (!linear) m = 1 + RND.irandom(8);
    boost::tuple<Polynomial<kind>,Polynomial<kind>,unsigned int> triple(alpha,beta,m);
    terms.push_back(triple);
  }
  if (!homogeneous) {
    Polynomial<kind> alpha;
    remainder = alpha;
  }
}

template<class kind>
Variety<unsigned int> Functional_Equation<kind>::reduce(unsigned int p)
{
  unsigned int i,j,in1;
  std::pair<unsigned int,unsigned int> duo;
  boost::tuple<Polynomial<kind>,Polynomial<kind>,unsigned int> trio;
  Variety<unsigned int> output(p,p);
  Polynomial<kind> py;
  Monomial<unsigned int> term;

  output.clear();

  for(i=0; i<p; ++i) {
    for(j=0; j<terms.size(); ++j) {
      trio = terms[j];
      // First convert alpha(p) to an element of GF(p)
      py = boost::get<0>(trio);
      in1 = convert(py.evaluate(i),p);
      if (in1 == 0) continue;
      term.coefficient = in1;
      // Next we need to convert beta(p) to an element of GF(p)
      py = boost::get<1>(trio);
      duo.first = convert(py.evaluate(i),p);
      // Finally we have to find out the exponent 
      duo.second = boost::get<2>(trio);
      term.exponents.push_back(duo);
      output.add_term(i,term);
      term.exponents.clear();
    }
    in1 = convert(remainder.evaluate(i),p);
    output.set_value(i,in1);
  }
  output.elaborate();
  return output;
}
