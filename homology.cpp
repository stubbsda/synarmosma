#include "homology.h"

using namespace SYNARMOSMA;

extern template const NTL::ZZ Integer_Matrix<NTL::ZZ>::neg1;
extern template const NTL::ZZ Integer_Matrix<NTL::ZZ>::unity; 

Homology::Homology()
{

}

Homology::Homology(Field f,Method m)
{
  field = f;
  method = m;
}

Homology::Homology(const Homology& source)
{
  betti_number = source.betti_number;
  torsion = source.torsion;
  field = source.field;
  method = source.method;
}

Homology& Homology::operator =(const Homology& source)
{
  if (this == &source) return *this;

  betti_number = source.betti_number;
  torsion = source.torsion;
  field = source.field;
  method = source.method;

  return *this;
}

Homology::~Homology()
{

}

void Homology::clear()
{
  field = Field::int32;
  method = Method::native;
  betti_number.clear();
  torsion.clear();
}

std::string Homology::write() const
{ 
  std::string output = "{";

  if (betti_number.empty()) {
    output += "}";
    return output;
  }
  unsigned int i,j,m,n = betti_number.size() - 1;
  std::stringstream sstream;
  sstream << "{";
  for(i=0; i<n; ++i) {
    sstream << "[" << betti_number[i];
    if (torsion[i].empty()) {
      sstream << "],";
      continue;
    } 
    sstream << ";";
    m = torsion[i].size() - 1;
    for(j=0; j<m; ++j) {
      sstream << torsion[i][j] << ",";
    }
    sstream << torsion[i][m] << "],";
  }
  sstream << "[" << betti_number[n];
  if (torsion[n].empty()) {
    sstream << "]}";
    return sstream.str();
  } 
  m = torsion[n].size() - 1;
  for(j=0; j<m; ++j) {
    sstream << torsion[n][j] << ",";
  }
  sstream << torsion[n][m] << "]}";  

  return sstream.str();
}

int Homology::serialize(std::ofstream& s) const 
{
  unsigned int i,j,k,m,n = betti_number.size();
  int count = 0;

  s.write((char*)(&field),sizeof(Field)); count += sizeof(Field);
  s.write((char*)(&method),sizeof(Method)); count += sizeof(Method);
  s.write((char*)(&n),sizeof(int)); count += sizeof(int);
  for(i=0; i<n; ++i) {
    m = betti_number[i];
    s.write((char*)(&m),sizeof(int)); count += sizeof(int);
    m = torsion[i].size();
    s.write((char*)(&m),sizeof(int)); count += sizeof(int);
    for(j=0; j<m; ++j) {
      k = torsion[i][j];
      s.write((char*)(&k),sizeof(int)); count += sizeof(int);
    }
  }
  return count;
}

int Homology::deserialize(std::ifstream& s)
{
  unsigned int i,j,k,n,m;
  int count = 0;
  std::vector<unsigned int> tau;

  clear();

  s.read((char*)(&field),sizeof(Field)); count += sizeof(Field);
  s.read((char*)(&method),sizeof(Method)); count += sizeof(Method);
  s.read((char*)(&n),sizeof(int)); count += sizeof(int);
  for(i=0; i<n; ++i) {
    s.read((char*)(&m),sizeof(int)); count += sizeof(int);
    betti_number.push_back(m);
    s.read((char*)(&m),sizeof(int)); count += sizeof(int);
    for(j=0; j<m; ++j) {
      s.read((char*)(&k),sizeof(int)); count += sizeof(int);
      tau.push_back(k);
    }
    torsion.push_back(tau);
    tau.clear();
  }
  return count;
} 

void Homology::compute_integral_native(const Nexus* NX)
{
#ifdef DEBUG
  assert(field == Field::int32 || field == Field::multiprecision);
#endif
  const int dimension = NX->get_dimension();
  const int nvertex = NX->get_order();
  int i,j,d,p,v2[2];
  unsigned int k,r,ulimit,d1 = nvertex,d2 = NX->get_length(1);
  std::string cx;
  std::stringstream ss;
  std::vector<unsigned int> image,kernel,tgenerators;
  std::set<int> vx,vy,faces,S;
  std::set<int>::const_iterator it;
  //hash_map::const_iterator qt;

  image.push_back(0);
  kernel.push_back(0);
  for(d=1; d<=dimension; ++d) {
    betti_number.push_back(0);
    torsion.push_back(tgenerators);
  } 

  if (field == Field::int32) {
    int alpha;
    Integer_Matrix<int>* A = new Integer_Matrix<int>(d1,d2);
    for(i=0; i<nvertex; ++i) {
      NX->get_neighbours(i,vx);
      for(it=vx.begin(); it!=vx.end(); ++it) {
        p = *it;
        S.insert(i); S.insert(p);
        j = NX->get_index(S);
        //qt = NX->index_table[1].find(S);
        //j = qt->second;
        S.clear();
        NX->get_elements(1,j,v2);
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
      d1 = NX->get_length(d-1);
      d2 = NX->get_length(d);
      A->initialize(d1,d2);
      for(k=0; k<d1; ++k) {
        NX->get_elements(d-1,k,vx);
        NX->get_entourage(d-1,k,faces);
        for(it=faces.begin(); it!=faces.end(); ++it) {
          j = *it;
          NX->get_elements(d,j,vy);
          alpha = coincidence(vx,vy);
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
    Integer_Matrix<NTL::ZZ>* A = new Integer_Matrix<NTL::ZZ>(d1,d2);

    for(i=0; i<nvertex; ++i) {
      NX->get_neighbours(i,vx);
      for(it=vx.begin(); it!=vx.end(); ++it) {
        p = *it;
        S.insert(i);
        S.insert(p);
        j = NX->get_index(S);
        NX->get_elements(1,j,v2);
        S.clear();
        if (i == v2[0]) {
          A->set(i,j,Integer_Matrix<NTL::ZZ>::neg1); 
        }
        else if (i == v2[1]) {
          A->set(i,j,Integer_Matrix<NTL::ZZ>::unity); 
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
      if (alpha > Integer_Matrix<NTL::ZZ>::unity) {
        ss << alpha;
        cx = ss.str();
        tgenerators.push_back(boost::lexical_cast<int>(cx.c_str()));
        ss.str("");
      }
    }
    torsion[1] = tgenerators;
    tgenerators.clear();    

    for(d=2; d<=dimension; ++d) {
      d1 = NX->get_length(d-1);
      d2 = NX->get_length(d);
      A->initialize(d1,d2);
      for(k=0; k<d1; ++k) {
        NX->get_elements(d-1,k,vx);
        NX->get_entourage(d-1,k,faces);
        for(it=faces.begin(); it!=faces.end(); ++it) {
          j = *it;
          NX->get_elements(d,j,vy);
          alpha = NTL::to_ZZ(coincidence(vx,vy));
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
        if (alpha > Integer_Matrix<NTL::ZZ>::unity) {
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
    betti_number[d] = kernel[d] - image[d+1];
  }
  betti_number[dimension] = kernel[dimension];
}

void Homology::compute_native(const Nexus* NX)
{
  int d;
  unsigned int betti = 1;
  std::vector<unsigned int> null;

  betti_number.clear();
  torsion.clear();

  if (!NX->connected()) {
    // In this case the integral homology group is just the free abelian group
    // on the number of distinct components...
    std::vector<int> components;
    betti = NX->component_analysis(components);
  }
  betti_number.push_back(betti);
  torsion.push_back(null);

  // Now the higher order homology groups...
  if (field == Field::mod2) {
    // Note that torsion isn't possible in this case - all we need to do is compute the Betti numbers
    int i,j,p,alpha;
    unsigned int r,k,d1 = NX->get_order(),d2 = NX->get_length(1); 
    std::vector<unsigned int> image,kernel;
    std::set<int> vx,vy,faces,S;
    std::set<int>::const_iterator it;
    Binary_Matrix* A = new Binary_Matrix(d1,d2);

    image.push_back(0);
    kernel.push_back(0);

    for(i=0; i<NX->get_order(); ++i) {
      NX->get_neighbours(i,vx);
      for(it=vx.begin(); it!=vx.end(); ++it) {
        p = *it;
        S.insert(i);
        S.insert(p);
        j = NX->get_index(S);
        S.clear();
        A->set(i,j);
      }
    }
    r = A->rank();
    image.push_back(r);
    kernel.push_back(d2 - r);

    for(d=2; d<=NX->get_dimension(); ++d) {
      d1 = NX->get_length(d-1);
      d2 = NX->get_length(d); 
      A->initialize(d1,d2);
      for(k=0; k<d1; ++k) {
        NX->get_elements(d-1,k,vx);
        NX->get_entourage(d-1,k,faces);
        for(it=faces.begin(); it!=faces.end(); ++it) {
          j = *it;
          NX->get_elements(d,j,vy);
          alpha = coincidence(vx,vy);
          if (alpha != 0) A->set(k,j);   
        }
      }
      r = A->rank();
      image.push_back(r);
      kernel.push_back(d2 - r);
    }
    for(d=1; d<NX->get_dimension(); ++d) {
      betti_number.push_back(kernel[d] - image[d+1]);
      torsion.push_back(null);
    }
    betti_number.push_back(kernel[NX->get_dimension()]);
    torsion.push_back(null);
    delete A;
  }
  else {
    compute_integral_native(NX);
  }
}

void Homology::compute_gap(const Nexus* NX)
{
  if (NX->get_dimension() < 1) return;
  int n;
  std::string line;
  std::set<int>::const_iterator it;

  betti_number.clear();
  torsion.clear();

  if (field == Field::mod2) {
    // Construct the boundary operator at each dimension and use GAP to 
    // compute the rank of this matrix over GF2...
    unsigned int i,j,d1,d2,r,betti = 1;
    bool done;
    int d,alpha;
    std::vector<unsigned int> null,image,kernel;
    std::vector<std::string> tokens;
    std::vector<unsigned char> rvector;
    std::set<int> vx,vy,faces;
    
    kernel.push_back(0);
    image.push_back(0);
    if (!(NX->connected())) {
      // In this case the integral homology group is just the free abelian group
      // on the number of distinct components...
      std::vector<int> components;
      betti = NX->component_analysis(components);
    }
    betti_number.push_back(betti);
    torsion.push_back(null);

    for(d=1; d<=NX->get_dimension(); ++d) {
      std::ofstream s("input.gap");
      s << "A := [";
      d1 = NX->get_length(d-1);
      d2 = NX->get_length(d);
      for(i=0; i<d1; ++i) {        
        for(j=0; j<d2; ++j) {
          rvector.push_back(0);
        }
        NX->get_elements(d-1,i,vx);
        NX->get_entourage(d-1,i,faces);
        for(it=faces.begin(); it!=faces.end(); ++it) {
          j = *it;
          NX->get_elements(d,j,vy);
          alpha = coincidence(vx,vy);
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
      if (std::system("gap -b < input.gap > homology.dat") != 0) throw std::runtime_error("Error in GAP execution!");
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
    for(d=1; d<NX->get_dimension(); ++d) {
      betti_number.push_back(kernel[d] - image[d+1]);
      torsion.push_back(null);
    }
    betti_number.push_back(kernel[NX->get_dimension()]);
    torsion.push_back(null);
    return;
  }

  int i,j,k;
  std::vector<unsigned int> tnumber;
  std::set<int> vx,faces;
  bool first = true;
  int depth;
  std::string hdata;
  std::string::size_type schar;

  std::ofstream s("input.gap");
  s << "LoadPackage(\"simpcomp\");;" << std::endl;
  s << "complex := SCFromFacets([";
  for(i=1; i<=NX->get_dimension(); ++i) {
    n = (signed) NX->get_length(i);
    k = 0;
    for(j=0; j<n; ++j) {
      NX->get_entourage(i,j,faces);
      if (!(faces.empty())) continue;
      if (first) {
        s << "[";
        first = false;
      }
      else {
        s << ",[";
      }
      k = 0;
      NX->get_elements(i,j,vx);
      for(it=vx.begin(); it!=vx.end(); ++it) {
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
  if (std::system("gap -b < input.gap > homology.dat") != 0) throw std::runtime_error("Error in GAP execution!");
  std::ifstream file("homology.dat",std::ios::in);
  while(std::getline(file,line)) {
    if (line.empty()) continue;
    schar = line.find("[");
    if (schar == std::string::npos) continue;
    hdata = line.substr(schar);
    break;      
  }
  file.close();
#ifdef DEBUG
  assert(!(hdata.empty()));
#endif
  // Now we need to parse the line that contains the homology groups' structure
  boost::char_separator<char> sep(", ");
  boost::tokenizer<boost::char_separator<char> > tok(hdata,sep);
  depth = 0;
  std::string str;
  bool analyze_torsion = false;
  bool group = false;

  for(boost::tokenizer<boost::char_separator<char> >::iterator beg=tok.begin(); beg!=tok.end(); beg++) {
    str = *beg;
    if (str == "[") depth++;
    if (str == "]") depth--;
    if (depth == 2 && !group) {
      group = true;
    }
    else if (depth == 2 && group && !analyze_torsion) {
      betti_number.push_back(boost::lexical_cast<int>(str));
    }
    else if (depth == 2) {
      analyze_torsion = false;
      torsion.push_back(tnumber);
      tnumber.clear();
    }
    if (depth == 3 && analyze_torsion) {
      tnumber.push_back(boost::lexical_cast<int>(str));
    }
    else if (depth == 3) {
      analyze_torsion = true;
    }
    if (depth == 1 && group) group = false;
  }
  // Due to a weird issue with simpcomp...
  betti_number[0] += 1;
}

void Homology::compute(const Nexus* NX)
{
  if (method == Method::gap) {
    compute_gap(NX);
  }
  else {
    compute_native(NX);
  }
}

