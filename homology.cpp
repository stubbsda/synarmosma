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

#include "homology.h"

using namespace SYNARMOSMA;

Homology::Homology()
{

}

Homology::Homology(FIELD f,METHOD m)
{
  field = f;
  method = m;
}

Homology::~Homology()
{

}

void Homology::betti_numbers(std::vector<int>& bnumbers) const
{
  unsigned int i;
  bnumbers.clear();
  for(i=0; i<sequence.size(); ++i) {
    bnumbers.push_back(sequence[i].get_rank());
  }
}

void Homology::append(const Group& g)
{
  sequence.push_back(g);
}

void Homology::clear()
{
  sequence.clear();
}

std::string Homology::write() const
{
  std::string output = "[";
  if (sequence.empty()) {
    output += "]";
    return output;
  }
  int i,n = (signed) sequence.size() - 1;
  
  for(i=0; i<n; ++i) {
    output += sequence[i].compact_form() + "],[";
  }
  output += sequence[n].compact_form() + "]";

  return output;
}

void Homology::serialize(std::ofstream& s) const 
{
  int i,n = (signed) sequence.size();

  s.write((char*)(&field),sizeof(FIELD));
  s.write((char*)(&method),sizeof(METHOD));
  s.write((char*)(&n),sizeof(int));
  for(i=0; i<n; ++i) {
    sequence[i].serialize(s);
  }
}

void Homology::deserialize(std::ifstream& s)
{
  int i,n;
  Group g;

  clear();

  s.read((char*)(&field),sizeof(FIELD));
  s.read((char*)(&method),sizeof(METHOD));
  s.read((char*)(&n),sizeof(int));
  for(i=0; i<n; ++i) {
    g.deserialize(s);
    sequence.push_back(g);
  }  
} 

void Homology::compute_integral_native(const Nexus* NX)
{
  assert(field == INT || field == ZZ);
  const int dimension = NX->dimension;
  const int nvertex = NX->nvertex;
  int i,j,d,p,v2[2];
  unsigned int k,r,betti,ulimit,d1 = nvertex,d2 = NX->elements[1].size();
  std::string cx;
  std::stringstream ss;
  std::vector<unsigned int> image,kernel,tgenerators;
  std::vector<unsigned int>* torsion = new std::vector<unsigned int>[1+dimension];
  std::set<int> vx,S;
  std::set<int>::const_iterator it;
  hash_map::const_iterator qt;

  image.push_back(0);
  kernel.push_back(0);
  for(d=0; d<=dimension; ++d) {
    torsion[d] = tgenerators;
  } 

  if (field == INT) {
    int alpha;
    Matrix<int>* A = new Matrix<int>(d1,d2);

    for(i=0; i<nvertex; ++i) {
      for(it=NX->neighbours[i].begin(); it!=NX->neighbours[i].end(); ++it) {
        p = *it;
        S.insert(i); S.insert(p);
        qt = NX->index_table[1].find(S);
        j = qt->second;
        S.clear();
        NX->elements[1][j].get_vertices(v2);
        if (i == v2[0]) {
          A->set(i,j,-1); 
        }
        else if (i == v2[1]) {
          A->set(i,j,1); 
        }
      }
    }
    normalize(*A);
    r = 0;
    for(k=0; k<d1; ++k) {
      if (A->empty_row(k)) continue;
      r++;
    }
    image.push_back(r);
    kernel.push_back(d2 - r);
    ulimit = (d1 < d2) ? d1 : d2;
    for(k=0; k<ulimit; ++k) {
      if (A->empty_row(k)) continue;
      alpha = A->get_first_nonzero(k);
      if (alpha > 1) {
        ss << alpha;
        cx = ss.str();
        tgenerators.push_back(boost::lexical_cast<int>(cx.c_str()));
        ss.str("");
      }
    }
    torsion[1] = tgenerators;
    tgenerators.clear();

    for(d=2; d<=dimension; ++d) {
      d1 = NX->elements[d-1].size();
      d2 = NX->elements[d].size();
      A->initialize(d1,d2);
      for(k=0; k<d1; ++k) {
        vx = NX->elements[d-1][k].vertices;
        for(it=NX->elements[d-1][k].entourage.begin(); it!=NX->elements[d-1][k].entourage.end(); ++it) {
          j = *it;
          alpha = coincidence(vx,NX->elements[d][j].vertices);
          if (alpha != 0) A->set(k,j,alpha); 
        }
      }
      normalize(*A);
      r = 0;
      for(k=0; k<d1; ++k) {
        if (A->empty_row(k)) continue;
        r++;
      }
      image.push_back(r);
      kernel.push_back(d2 - r);
      ulimit = (d1 < d2) ? d1 : d2;
      for(k=0; k<ulimit; ++k) {
        if (A->empty_row(k)) continue;
        alpha = A->get_first_nonzero(k);
        if (alpha > 1) {
          ss << alpha;
          cx = ss.str();
          tgenerators.push_back(boost::lexical_cast<int>(cx.c_str()));
          ss.str("");
        }
      }
      torsion[d] = tgenerators;
      tgenerators.clear();
    }
    delete A;
  }
  else {
    NTL::ZZ alpha;
    Matrix<NTL::ZZ>* A = new Matrix<NTL::ZZ>(d1,d2);

    for(i=0; i<nvertex; ++i) {
      for(it=NX->neighbours[i].begin(); it!=NX->neighbours[i].end(); ++it) {
        p = *it;
        S.insert(i);
        S.insert(p);
        qt = NX->index_table[1].find(S);
        j = qt->second;
        NX->elements[1][j].get_vertices(v2);
        S.clear();
        if (i == v2[0]) {
          A->set(i,j,Matrix<NTL::ZZ>::neg1); 
        }
        else if (i == v2[1]) {
          A->set(i,j,Matrix<NTL::ZZ>::unity); 
        }
      }
    }
    normalize(*A);
    r = 0;
    for(k=0; k<d1; ++k) {
      if (A->empty_row(k)) continue;
      r++;
    }
    image.push_back(r);
    kernel.push_back(d2 - r);
    ulimit = (d1 < d2) ? d1 : d2;
    for(k=0; k<ulimit; ++k) {
      if (A->empty_row(k)) continue;
      alpha = A->get_first_nonzero(k);
      if (alpha > Matrix<NTL::ZZ>::unity) {
        ss << alpha;
        cx = ss.str();
        tgenerators.push_back(boost::lexical_cast<int>(cx.c_str()));
        ss.str("");
      }
    }
    torsion[1] = tgenerators;
    tgenerators.clear();    

    for(d=2; d<=dimension; ++d) {
      d1 = NX->elements[d-1].size();
      d2 = NX->elements[d].size();
      A->initialize(d1,d2);
      for(k=0; k<d1; ++k) {
        vx = NX->elements[d-1][k].vertices;
        for(it=NX->elements[d-1][k].entourage.begin(); it!=NX->elements[d-1][k].entourage.end(); ++it) {
          j = *it;
          alpha = NTL::to_ZZ(coincidence(vx,NX->elements[d][j].vertices));
          if (alpha != 0) A->set(k,j,alpha); 
        }
      }
      normalize(*A);
      r = 0;
      for(k=0; k<d1; ++k) {
        if (A->empty_row(k)) continue;
        r++;
      }
      image.push_back(r);
      kernel.push_back(d2 - r);
      ulimit = (d1 < d2) ? d1 : d2;
      for(k=0; k<ulimit; ++k) {
        if (A->empty_row(k)) continue;
        alpha = A->get_first_nonzero(k);
        if (alpha > Matrix<NTL::ZZ>::unity) {
          ss << alpha;
          cx = ss.str();
          tgenerators.push_back(boost::lexical_cast<int>(cx.c_str()));
          ss.str("");
        }
      }
      torsion[d] = tgenerators;
      tgenerators.clear();
    }
    delete A;
  }
  for(d=1; d<dimension; ++d) {
    betti = kernel[d] - image[d+1];
    sequence.push_back(Group(betti,torsion[d])); 
  }
  betti = kernel[dimension];
  sequence.push_back(Group(betti,torsion[dimension]));
  delete[] torsion;
}

void Homology::compute_native(const Nexus* NX)
{
  int d;
  unsigned int betti = 1;
  std::vector<unsigned int> torsion;

  sequence.clear();

  if (!NX->connected()) {
    // In this case the integral homology group is just the free abelian group
    // on the number of distinct components...
    std::vector<int> components;
    betti = NX->component_analysis(components);
  }
  sequence.push_back(Group(betti,torsion));

  // Now the higher order homology groups...
  if (field == GF2) {
    // Note that torsion isn't possible in this case - all we need to do is compute the Betti numbers
    int i,j,p,alpha;
    unsigned int r,k,d1 = NX->nvertex,d2 = NX->elements[1].size();
    std::vector<unsigned int> image,kernel;
    std::set<int> vx,S;
    std::set<int>::const_iterator it;
    hash_map::const_iterator qt;
    Binary_Matrix* A = new Binary_Matrix(d1,d2);

    image.push_back(0);
    kernel.push_back(0);

    for(i=0; i<NX->nvertex; ++i) {
      for(it=NX->neighbours[i].begin(); it!=NX->neighbours[i].end(); ++it) {
        p = *it;
        S.insert(i);
        S.insert(p);
        qt = NX->index_table[1].find(S);
        j = qt->second;
        S.clear();
        A->set(i,j);
      }
    }
    r = A->rank();
    image.push_back(r);
    kernel.push_back(d2 - r);

    for(d=2; d<=NX->dimension; ++d) {
      d1 = NX->elements[d-1].size();
      d2 = NX->elements[d].size();
      A->initialize(d1,d2);
      for(k=0; k<d1; ++k) {
        vx = NX->elements[d-1][k].vertices;
        for(it=NX->elements[d-1][k].entourage.begin(); it!=NX->elements[d-1][k].entourage.end(); ++it) {
          j = *it;
          alpha = coincidence(vx,NX->elements[d][j].vertices);
          if (alpha != 0) A->set(k,j);   
        }
      }
      r = A->rank();
      image.push_back(r);
      kernel.push_back(d2 - r);
    }
    for(d=1; d<NX->dimension; ++d) {
      betti = kernel[d] - image[d+1];
      sequence.push_back(Group(betti,torsion));
    }
    betti = kernel[NX->dimension];
    sequence.push_back(Group(betti,torsion));
    delete A;
  }
  else {
    compute_integral_native(NX);
  }
}

void Homology::compute_gap(const Nexus* NX)
{
  if (NX->dimension < 1) return;
  int n;
  std::string line;
  std::set<int>::const_iterator it;

  sequence.clear();

  if (field == GF2) {
    // Construct the boundary operator at each dimension and use GAP to 
    // compute the rank of this matrix over GF2...
    unsigned int i,j,d1,d2,r,betti = 1;
    bool done;
    int d,alpha;
    std::vector<unsigned int> torsion,image,kernel;
    std::vector<std::string> tokens;
    std::vector<unsigned char> rvector;
    std::set<int> vx;
    
    kernel.push_back(0);
    image.push_back(0);
    if (!(NX->connected())) {
      // In this case the integral homology group is just the free abelian group
      // on the number of distinct components...
      std::vector<int> components;
      betti = NX->component_analysis(components);
    }
    sequence.push_back(Group(betti,torsion));

    for(d=1; d<=NX->dimension; ++d) {
      std::ofstream s("input.gap");
      s << "A := [";
      d1 = NX->elements[d-1].size();
      d2 = NX->elements[d].size();
      for(i=0; i<d1; ++i) {        
        for(j=0; j<d2; ++j) {
          rvector.push_back(0);
        }
        vx = NX->elements[d-1][i].vertices;
        for(it=NX->elements[d-1][i].entourage.begin(); it!=NX->elements[d-1][i].entourage.end(); ++it) {
          j = *it;
          alpha = coincidence(vx,NX->elements[d][j].vertices);
          if (alpha != 0) rvector[j] = 1;
        }
        s << "[";
        for(j=0; j<d2-1; ++j) {
          if (rvector[j] == 1) {
            s << "Z(2)^0,";
          }
          else {
            s << "0*Z(2),";
          }
        }
        if (rvector[d2-1] == 1) {
          s << "Z(2)^0";
        }
        else {
          s << "0*Z(2)";
        }
        if (i == d1-1) {
          s << "]";
        }
        else {
          s << "],";
        }
        rvector.clear();
      }
      s << "];;" << std::endl;
      s << "RankMat(A);" << std::endl;
      s << "quit;" << std::endl;
      s.close();
      if (std::system("gap -b < input.gap > homology.dat") != 0) {
        std::cerr << "Unable to execute the GAP software in Nexus::compute_homology, exiting..." << std::endl;
        std::exit(1);
      }
      // Now get the value of the matrix rank...
      std::ifstream file("homology.dat");
      while(std::getline(file,line)) {
        if (line.empty()) continue;
        // Skip if it's just a line identifying GAP and its version number etc.
        if (line.find("GAP") != std::string::npos) continue;
        // Break up the line on spaces and look for a number, this will be rank
        boost::split(tokens,line,boost::is_any_of(" "));
        done = false;
        for(i=0; i<tokens.size(); ++i) {
          if (tokens[i].empty()) continue;
          try {
            r = boost::lexical_cast<int>(tokens[i]);
            done = true;
          }
          catch(boost::bad_lexical_cast&) {
            continue;
          }
          if (done) break;
        }
      }
      file.close();
      image.push_back(r);
      kernel.push_back(d2 - r);
    }
    for(d=1; d<NX->dimension; ++d) {
      betti = kernel[d] - image[d+1];
      sequence.push_back(Group(betti,torsion));
    }
    betti = kernel[NX->dimension];
    sequence.push_back(Group(betti,torsion));
    return;
  }

  int i,j,k;
  std::vector<Word> relations;
  Word w(0);
  std::vector<unsigned int> bnumber,tnumber;
  std::vector<std::vector<unsigned int> > telements;
  bool first = true;
  int depth;
  std::string hdata;
  std::string::size_type schar;

  std::ofstream s("input.gap");
  s << "LoadPackage(\"simpcomp\");;" << std::endl;
  s << "complex := SCFromFacets([";
  for(i=1; i<=NX->dimension; ++i) {
    n = (signed) NX->elements[i].size();
    k = 0;
    for(j=0; j<n; ++j) {
      if (!(NX->elements[i][j].entourage.empty())) continue;
      if (first) {
        s << "[";
        first = false;
      }
      else {
        s << ",[";
      }
      k = 0;
      for(it=NX->elements[i][j].vertices.begin(); it!=NX->elements[i][j].vertices.end(); ++it) {
        if (k < i) {
          s << *it << ",";
        }
        else {
          s << *it;
        }
        ++k;
      }
      s << "]";
    }
  }
  s << "]);;";
  s << std::endl;
  s << "SCHomology(complex);" << std::endl;
  s << "quit;" << std::endl;
  s.close();
  if (std::system("gap -b < input.gap > homology.dat") != 0) {
    std::cerr << "Error in executing the GAP software in Nexus::compute_homology - are GAP and its simpcomp package installed?" << std::endl;
    std::cerr << "Exiting..." << std::endl;
    std::exit(1);
  }
  std::ifstream file("homology.dat",std::ios_base::in);
  while(std::getline(file,line)) {
    if (line.empty()) continue;
    schar = line.find("[");
    if (schar == std::string::npos) continue;
    hdata = line.substr(schar);
    break;      
  }
  file.close();
  assert(!(hdata.empty()));
  // Now we need to parse the line that contains the homology groups' structure
  boost::char_separator<char> sep(", ");
  boost::tokenizer<boost::char_separator<char> > tok(hdata,sep);
  depth = 0;
  std::string str;
  bool torsion = false;
  bool group = false;

  for(boost::tokenizer<boost::char_separator<char> >::iterator beg=tok.begin(); beg!=tok.end(); beg++) {
    str = *beg;
    if (str == "[") depth++;
    if (str == "]") depth--;
    if (depth == 2 && !group) {
      group = true;
    }
    else if (depth == 2 && group && !torsion) {
      bnumber.push_back(boost::lexical_cast<int>(str));
    }
    else if (depth == 2) {
      torsion = false;
      telements.push_back(tnumber);
      tnumber.clear();
    }
    if (depth == 3 && torsion) {
      tnumber.push_back(boost::lexical_cast<int>(str));
    }
    else if (depth == 3) {
      torsion = true;
    }
    if (depth == 1 && group) group = false;
  }
  // Due to a weird issue with simpcomp...
  bnumber[0] += 1;
  for(i=0; i<=NX->dimension; ++i) {
    sequence.push_back(Group(bnumber[i],telements[i]));
  }
}

void Homology::compute(const Nexus* NX)
{
  if (method == GAP) {
    compute_gap(NX);
  }
  else {
    compute_native(NX);
  }
}

