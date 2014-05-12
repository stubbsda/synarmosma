#include "propositional_system.h"

Propositional_System::Propositional_System()
{
  set_default_values();
  initialize(100);
}

Propositional_System::Propositional_System(unsigned int n)
{
  set_default_values();
  initialize(n);
}

Propositional_System::Propositional_System(unsigned int n,const char* filename)
{
  unsigned int i,eq_point;
  std::string line,name,value;

  set_default_values();

  // Open the parameter file
  std::ifstream s(filename,std::ios_base::in);
  if (s.is_open() == false) {
    // If the file isn't there, print an error message and
    // exit cleanly
    std::cout << "The file " << filename << " cannot be found!" << std::endl;
    std::exit(1);
  }

  while(std::getline(s,line)) {
    // If it's an empty line, continue
    if (line.empty()) continue;
    // If the line begins with a #, ignore it
    if (line[0] == '#') continue;
    // Find the position of the equals sign
    eq_point = 0;
    for(i=0; i<line.length(); ++i) {
      if (line[i] == '=') {
        eq_point = i;
        break;
      }
    }
    // If there's no equals sign in this line, continue
    if (eq_point < 1) continue;
    name = line.substr(0,eq_point-1);
    trim(name);
    value = line.substr(eq_point+1,line.length());
    trim(value);
    // Now that we have the parameter name, see if it matches
    // any of the known parameters. If so, read in the value and
    // assign it
    if (name == "naxiom") {
      natom = boost::lexical_cast<int>(value);
    }
  }
  initialize(n);
}

Propositional_System::Propositional_System(unsigned int n,unsigned int m)
{
  // n = axiomezahl und m = natom
  set_default_values();
  natom = m;
  initialize(n);
}

Propositional_System::~Propositional_System()
{

}

void Propositional_System::set_default_values()
{
  natom = 10;
}

void Propositional_System::initialize(unsigned int n)
{
  unsigned int i;

  for(i=0; i<n; ++i) {
    theorems.push_back(Proposition(natom));
  }
  nuniverse = ipow(2,natom);
  for(i=0; i<theorems.size(); ++i) {
    truth.push_back(boost::dynamic_bitset<>(nuniverse));
  }
  build_logical_universe();
  std::cout << "Computing the internal logic..." << std::endl;
  compute_internal_logic();
  for(i=0; i<nuniverse; ++i) {
    delete[] logical_universe[i];
  }
  delete[] logical_universe;
}

void Propositional_System::write(const char* filename,unsigned int pointer)
{
  int i;
  std::ofstream s(filename,std::ios_base::out | std::ios_base::ate | std::ios_base::binary);
  s.seekp(pointer);
  for(i=0; i<theorems.size(); ++i) {
    theorems[i].serialize(s);
  }
  s.close();
}

void Propositional_System::write(const char* filename)
{
  unsigned int i;
  std::ofstream s(filename,std::ios_base::out | std::ios_base::binary);
  s.seekp(0);
  for(i=0; i<theorems.size(); ++i) {
    theorems[i].serialize(s);
  }
  s.close();
}

void Propositional_System::read(const char* filename,unsigned int pointer)
{
  unsigned int i;
  std::ifstream s(filename,std::ios_base::in | std::ios_base::binary);
  s.seekg(pointer);
  for(i=0; i<theorems.size(); ++i) {
    theorems[i].deserialize(s);
  }
  // Close the file
  s.close();
}

void Propositional_System::read(const char* filename)
{
  unsigned int i;
  std::ifstream s(filename,std::ios_base::in | std::ios_base::binary);
  s.seekg(0);
  for(i=0; i<theorems.size(); ++i) {
    theorems[i].deserialize(s);
  }
  // Close the file
  s.close();
}

void Propositional_System::build_logical_universe()
{
  // A method that carries out the recursion once and for all for
  // creating a large matrix of all possible logical values for a
  // system of natomic propositions
  unsigned int i,j,p,idiv,in1;

  logical_universe = new bool*[nuniverse];
  for(i=0; i<nuniverse; ++i) {
    logical_universe[i] = new bool[natom];
  }

  for(i=0; i<nuniverse; ++i) {
    in1 = i;
    for(j=0; j<natom; ++j) {
      p = ipow(2,natom-(j-1));
      idiv = in1/p;
      logical_universe[i][j] = (idiv % 2) ? false: true;
      in1 -= p*idiv;
    }
  }
}

unsigned int Propositional_System::bit_count(unsigned int i) const
{
  return truth[i].count();
}

unsigned int Propositional_System::consistency(unsigned int i,unsigned int j,const std::string& op) const
{
  boost::dynamic_bitset<> temp(nuniverse);
  if (op == "and") {
    temp = truth[i] & truth[j];
  }
  else if (op == "xor") {
    temp = truth[i] ^ truth[j];
  }
  else {
    temp = truth[i] | truth[j];
  }
  return temp.count();
}

bool Propositional_System::implication(const Proposition& test) const
{
  return true;
}

bool Propositional_System::implication(unsigned int in1,const std::vector<unsigned int>& axioms) const
{
  unsigned int i;
  boost::dynamic_bitset<> temp(nuniverse);

  temp = truth[in1];
  for(i=0; i<axioms.size(); ++i) {
    temp |= ~truth[axioms[i]];
  }
  if (temp.count() == nuniverse) return true;
  return false;
}

void Propositional_System::compute_internal_logic()
{
  // The same basic idea as above, fill up the array of edge weights
  // with values but based on propositional logic rather than number
  // theory
  unsigned int i,j,n = theorems.size();
  bool out;

  // Loop through all vertices, if the topology has been modified, then
  // compute its truth value (do the set of propositions of its neighbours,
  // taken as axioms, imply this vertex's proposition as a logical
  // conclusion), which will be used to compute the edge weights

  for(i=0; i<n; ++i) {
    truth[i].reset();
    for(j=0; j<nuniverse; ++j) {
      out = theorems[i].evaluate(logical_universe[j]);
      if (out) truth[i].set(j);
    }
  }
}

void Propositional_System::compute_pairs(unsigned int* edges)
{
  // This method must determine which propositions imply other propositions among
  // our collection, "theorems".
  unsigned int i,j,kt;
  const unsigned int nnode = theorems.size();
  boost::dynamic_bitset<> temp(nuniverse),p1(nuniverse);
  unsigned int* sink = new unsigned int[nnode];
  unsigned int* source = new unsigned int[nnode];

  std::cout << "Computing the logical pairs..." << std::endl;

  const unsigned int nedge = (nnode*(nnode-1))/2;

  for(i=0; i<nedge; ++i) {
    edges[i] = 0;
  }
  for(i=0; i<nnode; ++i) {
    sink[i] = 1;
    source[i] = 1;
  }
  for(i=0; i<nnode; ++i) {
    p1 = truth[i];
    for(j=0; j<nnode; ++j) {
      if (i == j) continue;
      // p(i) implies p(j)? That means that (!p(j) v p(i)) is always true
      temp = ~truth[j] | p1;
      if (temp.count() == nuniverse) {
        // How to compute the edge index?
        if (j > i) {
          kt = (j-i-1) + nedge - ((nnode-i)*(nnode-i-1))/2;
          edges[kt] = (edges[kt] == 0) ? 1 : 3;
        }
        else {
          // We have already seen this edge!
          kt = (i-j-1) + nedge - ((nnode-j)*(nnode-j-1))/2;
          edges[kt] = (edges[kt] == 0) ? 2 : 3;
        }
        source[j] = 0;
        sink[i] = 0;
      }
    }
  }

  unsigned int null = 0;
  unsigned int ndirected = 0;
  unsigned int ndouble = 0;
  for(i=0; i<nedge; ++i) {
    if (edges[i] == 0) null++;
    if (edges[i] == 1 || edges[i] == 2) ndirected++;
    if (edges[i] == 3) ndouble++;
  }
  std::cout << "The percent of null edges is " << 100.0*double(null)/double(nedge) << std::endl;
  std::cout << "The percent of directed edges is " << 100.0*double(ndirected)/double(nedge) << std::endl;
  std::cout << "The percent of double edges is " << 100.0*double(ndouble)/double(nedge) << std::endl;
  // What is a logical sink? A proposition which doesn't imply any other propositions
  unsigned int nsink = 0;
  for(i=0; i<nnode; ++i) {
    // A proposition that never implies any other but is implied by at least one other...
    if (sink[i] == 1 && source[i] == 0) nsink++;
  }

  // What is a logical source? A proposition is never implied by other propositions
  unsigned int nsource = 0;
  for(i=0; i<nnode; ++i) {
    // A proposition that is never implied by any other but implies at least one other...
    if (source[i] == 1 && sink[i] == 0) nsource++;
  }
  if (nsource == 0 && nsink == 0) {
    std::cout << "There are no logical sources and sinks in this graph." << std::endl;
  }
  else if (nsource == 0) {
    std::cout << "There are no logical sources and " << nsink << " logical sinks in this graph." << std::endl;
  }
  else if (nsink == 0) {
    std::cout << "There are " << nsource << " logical sources and no logical sinks in this graph." << std::endl;
  }
  else {
    std::cout << "There are " << nsource << " logical sources and " << nsink << " logical sinks in this graph." << std::endl;
  }
}



