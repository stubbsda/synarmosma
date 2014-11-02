/*
  Copyright 2014 Daniel Stubbs

  This file is part of Synarmosma.

  Synarmosma is free software: you can redistribute it and/or modify 
  it under the terms of the GNU General Public License as published by 
  the Free Software Foundation, either version 3 of the License, or 
  (at your option) any later version.

  Synarmosma is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with Synarmosma.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "propositional_system.h"

using namespace SYNARMOSMA;

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
    std::cerr << "The file " << filename << " cannot be found!" << std::endl;
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
  // n = number of axioms und m = number of atoms
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
  compute_internal_logic();
  for(i=0; i<nuniverse; ++i) {
    delete[] logical_universe[i];
  }
  delete[] logical_universe;
}

void Propositional_System::write(const char* filename,unsigned int pointer)
{
  unsigned int i;
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

void Propositional_System::compute_pairs(unsigned int* edges,std::set<unsigned int>& sink,std::set<unsigned int>& source) const 
{
  // This method must determine which propositions imply other propositions among
  // our collection, "theorems".
  // The output scheme for the edges is as follows:
  // edges[i] = 0 => null edge
  // edges[i] = 1 or 2 => directed edge
  // edges[i] = 3 => multi-edge
  unsigned int i,j,kt;
  std::vector<unsigned int> vsink,vsource;
  boost::dynamic_bitset<> temp(nuniverse),p1(nuniverse);
  const unsigned int nnode = theorems.size();
  const unsigned int nedge = (nnode*(nnode-1))/2;

  for(i=0; i<nedge; ++i) {
    edges[i] = 0;
  }
  for(i=0; i<nnode; ++i) {
    vsink.push_back(1);
    vsource.push_back(1);
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
        vsource[j] = 0;
        vsink[i] = 0;
      }
    }
  }

  // What is a logical sink? A proposition which doesn't imply any other propositions
  sink.clear();
  for(i=0; i<nnode; ++i) {
    // A proposition that never implies any other but is implied by at least one other...
    if (vsink[i] == 1 && vsource[i] == 0) sink.insert(i);
  }

  // What is a logical source? A proposition is never implied by other propositions
  source.clear();
  for(i=0; i<nnode; ++i) {
    // A proposition that is never implied by any other but implies at least one other...
    if (vsource[i] == 1 && vsink[i] == 0) source.insert(i);
  }
}



