#include "nexus.h"

Nexus::Nexus() : Schema()
{
  dimension = -1;
}

Nexus::Nexus(int n) : Schema()
{
  elements = new std::vector<Cell>[1+n];
  index_table = new hash_map[1+n];
  dimension = n;
}

Nexus::~Nexus()
{
  if (dimension > -1) {
    delete[] elements;
    delete[] index_table;
  }
}

void Nexus::initialize(int d)
{
  clear();
  assert(d > 0);
  elements = new std::vector<Cell>[1+d];
  index_table = new hash_map[1+d];
  dimension = d;
  nvertex = 0;
}

void Nexus::initialize(int d,int nv)
{
  clear();
  assert(d > 0);
  assert(nv > 0);
  elements = new std::vector<Cell>[1+d];
  index_table = new hash_map[1+d];
  dimension = d;
  nvertex = nv;
  std::set<int> null;
  for(int i=0; i<nvertex; ++i) {
    neighbours.push_back(null);
  }
}

void Nexus::clear()
{
  if (dimension > -1) {
    delete[] elements;
    delete[] index_table;
  }
  neighbours.clear();
  dimension = -1;
  nvertex = 0;
}

void Nexus::paste(const std::set<int>& vx)
{
  Cell c(vx);
  int m = c.dimension();
  elements[m].push_back(c);
  index_table[m][c.key] = elements[m].size() - 1;
}

int Nexus::size() const
{
  return nvertex;
}

void Nexus::compute_neighbours()
{
  int i,vx[2];
  const int ne = (signed) elements[1].size();
  for(i=0; i<nvertex; ++i) {
    neighbours[i].clear();
  }
  for(i=0; i<ne; ++i) {
    elements[1][i].get_vertices(vx);
    neighbours[vx[1]].insert(vx[0]);
    neighbours[vx[0]].insert(vx[1]);
  }
}

void Nexus::compute_entourages()
{
  int i,j,k,ns;
  std::set<int> s;
  std::set<int>::const_iterator it;
  std::string fx;
  hash_map::const_iterator qt;

  for(i=dimension; i>=1; i--) {
    ns = (signed) elements[i].size();
    for(j=0; j<ns; ++j) {
      for(k=0; k<1+i; ++k) {
        fx = elements[i][j].faces[k];
        qt = index_table[i-1].find(fx);
        if (qt == index_table[i-1].end()) {
          std::cerr << "Entourage error: " << i << "  " << j << "  " << elements[i][j].key << "  "<< fx << std::endl;
          std::exit(1);
        }
        elements[i-1][qt->second].entourage.insert(j);
      }
    }
  }
}

void Nexus::regularization()
{
  // If a given d-simplex exists, we need to make sure that all of
  // its sub-simplices also exist in the nexus
  int i,j,k,n,m;
  Cell S;
  hash_map::const_iterator qt;

  for(i=dimension; i>=1; i--) {
    n = (signed) elements[i].size();
    m = (signed) elements[i-1].size();
    for(j=0; j<n; ++j) {
      for(k=0; k<1+i; ++k) {
        qt = index_table[i-1].find(elements[i][j].faces[k]);
        if (qt == index_table[i-1].end()) {
          S = Cell(elements[i][j].faces[k]);
          elements[i-1].push_back(S);
          index_table[i-1][S.key] = m;
          m++;
        }
      }
    }
  }
  compute_entourages();
  compute_neighbours();
}

void Nexus::closure(const std::set<std::string>& S,Nexus* NX,int* offset) const
{
  int i,n,m;
  unsigned int j;
  std::string key;
  std::set<int> vx;
  std::vector<Cell>* cell_array = new std::vector<Cell>[1+dimension];
  hash_map::const_iterator qt;
  std::set<std::string>::const_iterator st;
  std::set<int>::const_iterator it;

  for(st=S.begin(); st!=S.end(); st++) {
    key = *st;
    n = std::count(key.begin(),key.end(),':');
    if (n > 0) {
      qt = index_table[n].find(key);
      m = qt->second;
      cell_array[n].push_back(elements[n][m]);
    }
    else {
      vx.insert(boost::lexical_cast<int>(key));
    }
  }

  for(i=0; i<nvertex; ++i) {
    offset[i] = -1;
  }
  // Here n will be the number of vertices in the closure of the new Nexus while m 
  // will be its dimension
  n = 0;
  m = 0;
  for(i=1; i<=dimension; ++i) {
    if (cell_array[i].empty()) continue;
    m = i;
    for(j=0; j<cell_array[i].size(); ++j) {
      for(it=cell_array[i][j].vertices.begin(); it!=cell_array[i][j].vertices.end(); it++) {
        if (offset[*it] == -1) {
          offset[*it] = n;
          n++;
        }
      }
    }
  }
  for(it=vx.begin(); it!=vx.end(); it++) {
    if (offset[*it] == -1) {
      offset[*it] = n;
      n++;
    }
  }
  vx.clear();

  NX->initialize(m,n);
  for(i=1; i<=m; ++i) {
    if (cell_array[i].empty()) continue;
    for(j=0; j<cell_array[i].size(); ++j) {
      for(it=cell_array[i][j].vertices.begin(); it!=cell_array[i][j].vertices.end(); it++) {
        vx.insert(offset[*it]);
      }
      NX->paste(vx);
      vx.clear();
    }
  }
  NX->regularization();
  delete[] cell_array;
}

void Nexus::link(const std::set<std::string>& S,std::vector<Cell>* cell_array) const
{
  int i,offset1[nvertex],offset2[nvertex];
  unsigned int j,k;
  bool found;
  std::string key;
  std::vector<Cell>* cell_array1 = new std::vector<Cell>[1+dimension];
  std::vector<Cell>* cell_array2 = new std::vector<Cell>[1+dimension];
  std::set<int> vx;
  std::set<std::string> string_set;
  std::set<int>::const_iterator it;
  Nexus* N1 = new Nexus;
  Nexus* N2 = new Nexus;

  // The link of S is defined to be closure(star(S)) - star(closure(S)), we begin 
  // with the first term
  star(S,cell_array1);
  // Now we need to convert cell_array to a set of strings...
  for(i=0; i<=dimension; ++i) {
    for(j=0; j<cell_array1[i].size(); ++j) {
      string_set.insert(cell_array1[i][j].key);
    }
  }  
  closure(string_set,N1,offset1);

  // Now the second term...
  closure(S,N2,offset2);
  string_set.clear();
  for(i=0; i<=N2->dimension; ++i) {
    for(j=0; j<N2->elements[i].size(); ++j) {
      for(it=N2->elements[i][j].vertices.begin(); it!=N2->elements[i][j].vertices.end(); ++it) {
        vx.insert(offset2[*it]);
      }
      string_set.insert(make_key(vx));
      vx.clear();
    }
  }
  star(string_set,cell_array2);

  // Finally, construct the output by taking every cell that exists in N1 and, if it doesn't also 
  // exist in cell_array2, adding it to cell_array
  for(i=0; i<=N1->dimension; ++i) {
    for(j=0; j<N1->elements[i].size(); ++j) {
      for(it=N1->elements[i][j].vertices.begin(); it!=N1->elements[i][j].vertices.end(); ++it) {
        vx.insert(offset1[*it]);
      }
      key = make_key(vx);
      // Does this element exist somewhere in cell_array2?
      found = false;
      for(k=0; k<cell_array2[i].size(); ++k) {
        if (key == cell_array2[i][k].key) {
          found = true;
          break;
        }
      }
      vx.clear();
      if (found) continue;
      cell_array[i].push_back(N1->elements[i][j]);
    }    
  }

  delete[] cell_array1;
  delete[] cell_array2;
  delete N1; delete N2;
}

void Nexus::star(const std::set<std::string>& S,std::vector<Cell>* output) const
{
  int n,m,d,nc;
  std::string key;
  hash_map ckey;
  std::vector<Cell> clist;
  hash_map::const_iterator qt;
  std::set<std::string>::const_iterator it;

  for(it=S.begin(); it!=S.end(); it++) {
    key = *it;
    n = std::count(key.begin(),key.end(),':');
    if (n > 0) {
      qt = index_table[n].find(key);
      m = qt->second;
    }
    else {
      m = boost::lexical_cast<int>(key);
    }
    ascend(n,m,clist);
  }
  for(n=0; n<=dimension; ++n) {
    output[n].clear();
  }
  nc = (signed) clist.size();
  d = std::count(clist[0].key.begin(),clist[0].key.end(),':');
  output[d].push_back(clist[0]);
  ckey[clist[0].key] = 0;
  m = 1;
  for(n=1; n<nc; ++n) {
    qt = ckey.find(clist[n].key);
    if (qt == ckey.end()) {
      d = std::count(clist[n].key.begin(),clist[n].key.end(),':');
      output[d].push_back(clist[n]);
      ckey[clist[n].key] = m;
      m++;
    }
  }
}

void Nexus::ascend(int d,int in1,std::vector<Cell>& out) const
{
  if (d == 1+ND) return;

  int i;
  std::set<int>::const_iterator it;
  std::set<int> S = elements[d][in1].entourage;

  for(it=S.begin(); it!=S.end(); it++) {
    i = *it;
    out.push_back(elements[d+1][i]);
    ascend(1+d,i,out);
  }
}

void Nexus::compute_integral_homology(std::vector<Group>& HZ,FIELD base) const
{
  assert(base == INT || base == ZZ);
  int i,j,d,p,v2[2];
  unsigned int k,r,betti,ulimit,d1 = nvertex,d2 = elements[1].size();
  std::string cx;
  std::stringstream ss;
  std::vector<unsigned int> image,kernel,tgenerators;
  std::vector<unsigned int>* torsion = new std::vector<unsigned int>[dimension+1];
  std::set<int> vx;
  std::set<int>::const_iterator it;
  hash_map::const_iterator qt;

  image.push_back(0);
  kernel.push_back(0);
  for(d=0; d<=dimension; ++d) {
    torsion[d] = tgenerators;
  } 

  if (base == INT) {
    int alpha;
    Matrix<int>* A = new Matrix<int>(d1,d2);

    for(i=0; i<nvertex; ++i) {
      for(it=neighbours[i].begin(); it!=neighbours[i].end(); it++) {
        p = *it;
        qt = index_table[1].find(make_key(i,p));
        j = qt->second;
        elements[1][j].get_vertices(v2);
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
      d1 = elements[d-1].size();
      d2 = elements[d].size();
      A->initialize(d1,d2);
      for(k=0; k<d1; ++k) {
        vx = elements[d-1][k].vertices;
        for(it=elements[d-1][k].entourage.begin(); it!=elements[d-1][k].entourage.end(); it++) {
          j = *it;
          alpha = coincidence(vx,elements[d][j].vertices);
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
      for(it=neighbours[i].begin(); it!=neighbours[i].end(); it++) {
        p = *it;
        qt = index_table[1].find(make_key(i,p));
        j = qt->second;
        elements[1][j].get_vertices(v2);
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
      d1 = elements[d-1].size();
      d2 = elements[d].size();
      A->initialize(d1,d2);
      for(k=0; k<d1; ++k) {
        vx = elements[d-1][k].vertices;
        for(it=elements[d-1][k].entourage.begin(); it!=elements[d-1][k].entourage.end(); it++) {
          j = *it;
          alpha = NTL::to_ZZ(coincidence(vx,elements[d][j].vertices));
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
    HZ.push_back(Group(betti,torsion[d])); 
  }
  betti = kernel[dimension];
  HZ.push_back(Group(betti,torsion[dimension]));
  delete[] torsion;
}

void Nexus::compute_homology_native(std::vector<Group>& HZ,FIELD base) const
{
  int d;
  unsigned int betti = 1;
  std::vector<unsigned int> torsion;

  HZ.clear();

  if (!connected()) {
    // In this case the integral homology group is just the free abelian group
    // on the number of distinct components...
    std::vector<int> components;
    betti = component_analysis(components);
  }
  HZ.push_back(Group(betti,torsion));

  // Now the higher order homology groups...
  if (base == GF2) {
    // Note that torsion isn't possible in this case - all we need to do is compute the Betti numbers
    int i,j,p,alpha;
    unsigned int r,k,d1 = nvertex,d2 = elements[1].size();
    std::vector<unsigned int> image,kernel;
    std::set<int> vx;
    std::set<int>::const_iterator it;
    hash_map::const_iterator qt;
    Binary_Matrix* A = new Binary_Matrix(d1,d2);

    image.push_back(0);
    kernel.push_back(0);

    for(i=0; i<nvertex; ++i) {
      for(it=neighbours[i].begin(); it!=neighbours[i].end(); it++) {
        p = *it;
        qt = index_table[1].find(make_key(i,p));
        j = qt->second;
        A->set(i,j);
      }
    }
    r = A->rank();
    image.push_back(r);
    kernel.push_back(d2 - r);

    for(d=2; d<=dimension; ++d) {
      d1 = elements[d-1].size();
      d2 = elements[d].size();
      A->initialize(d1,d2);
      for(k=0; k<d1; ++k) {
        vx = elements[d-1][k].vertices;
        for(it=elements[d-1][k].entourage.begin(); it!=elements[d-1][k].entourage.end(); it++) {
          j = *it;
          alpha = coincidence(vx,elements[d][j].vertices);
          if (alpha != 0) A->set(k,j);   
        }
      }
      r = A->rank();
      image.push_back(r);
      kernel.push_back(d2 - r);
    }
    for(d=1; d<dimension; ++d) {
      betti = kernel[d] - image[d+1];
      HZ.push_back(Group(betti,torsion));
    }
    betti = kernel[dimension];
    HZ.push_back(Group(betti,torsion));
    delete A;
  }
  else {
    compute_integral_homology(HZ,base);
  }
}

void Nexus::compute_homology(std::vector<Group>& HZ,FIELD base) const
{
  if (dimension < 1) return;
  int n;
  std::string line;
  std::set<int>::const_iterator it;

  HZ.clear();

  if (base == GF2) {
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
    if (!connected()) {
      // In this case the integral homology group is just the free abelian group
      // on the number of distinct components...
      std::vector<int> components;
      betti = component_analysis(components);
    }
    HZ.push_back(Group(betti,torsion));

    for(d=1; d<=dimension; ++d) {
      std::ofstream s("input.gap");
      s << "A := [";
      d1 = elements[d-1].size();
      d2 = elements[d].size();
      for(i=0; i<d1; ++i) {        
        for(j=0; j<d2; ++j) {
          rvector.push_back(0);
        }
        vx = elements[d-1][i].vertices;
        for(it=elements[d-1][i].entourage.begin(); it!=elements[d-1][i].entourage.end(); it++) {
          j = *it;
          alpha = coincidence(vx,elements[d][j].vertices);
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
    for(d=1; d<dimension; ++d) {
      betti = kernel[d] - image[d+1];
      HZ.push_back(Group(betti,torsion));
    }
    betti = kernel[dimension];
    HZ.push_back(Group(betti,torsion));
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
  for(i=1; i<=dimension; ++i) {
    n = (signed) elements[i].size();
    k = 0;
    for(j=0; j<n; ++j) {
      if (!elements[i][j].entourage.empty()) continue;
      if (first) {
        s << "[";
        first = false;
      }
      else {
        s << ",[";
      }
      k = 0;
      for(it=elements[i][j].vertices.begin(); it!=elements[i][j].vertices.end(); it++) {
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
  for(i=0; i<=dimension; ++i) {
    HZ.push_back(Group(bnumber[i],telements[i]));
  }
}

bool Nexus::pseudomanifold(bool* boundary) const
{
  int i,j;
  bool found,output = true;
  std::vector<int> candidates;

  if (dimension == 0) return false;

  *boundary = false;

  // First the non-branching requirement...
  for(i=0; i<(signed) elements[dimension-1].size(); ++i) {
    j = (signed) elements[dimension-1][i].entourage.size();
    // If the entourage size is ever 1, this is a pseudomanifold with boundary!
    if (j > 2 || j < 1) return false;
    if (j == 1) *boundary = true;
  }
  // Next, dimensional homogeneity...
  for(i=0; i<nvertex; ++i) {
    found = false;
    for(j=0; j<(signed) elements[dimension].size(); ++j) {
      if (elements[dimension][j].contains(i)) {
        found = true;
        break;
      }
    }
    if (!found) return false;
  }
  // Finally, as this is a non-branching simplicial complex, next
  // check if it is strongly connected.
  for(i=0; i<(signed) elements[dimension].size(); ++i) {
    candidates.push_back(i);
  }

  if (candidates.size() == 2) {
    // The simplest case: we just need to verify that these
    // two i-simplices differ by a single vertex
    if (affinity(elements[dimension][candidates[0]],elements[dimension][candidates[1]]) != dimension) return false;
  }
  else if (candidates.size() > 2) {
    Graph* G = new Graph;
    int in1,k,l,nf,ni,off1,off2,vx[dimension+1];
    std::vector<int> vlist;
    const int nc = (signed) candidates.size();
    const int N = nc*(nc-1)/2;
    int* link = new int[N];

    for(i=0; i<nc; ++i) {
      G->add_vertex();
      elements[dimension][candidates[i]].get_vertices(vx);
      for(j=0; j<=dimension; ++j) {
        vlist.push_back(vx[j]);
      }
    }
#ifdef PARALLEL
#pragma omp parallel for default(shared) private(i,j,k,l,in1,ni,nf,off1,off2) schedule(dynamic,1)
#endif
    for(i=0; i<nc; ++i) {
      ni = i*nc - i*(i+1)/2;
      off1 = (1+dimension)*i;
      for(j=1+i; j<nc; ++j) {
        nf = 0;
        off2 = (1+dimension)*j;
        for(k=0; k<=dimension; ++k) {
          in1 = vlist[off1+k];
          for(l=0; l<=dimension; ++l) {
            if (in1 == vlist[off2+l]) {
              nf++;
              break;
            }
          }
        }
        link[ni+j-(i+1)] = (nf == dimension) ? 1 : 0;
      }
    }
    for(i=0; i<nc; ++i) {
      ni = i*nc - i*(i+1)/2;
      for(j=1+i; j<nc; ++j) {
        if (link[ni+j-(1+i)] == 1) G->add_edge(i,j);
      }
    }
    if (!G->connected()) output = false;
    G->clear();
    delete G;
    delete[] link;
  }
  return output;
}

bool Nexus::orientable() const
{
  // Program to calculate whether the complex is orientable, note
  // this only makes sense if the complex is also a pseudomanifold!

  // If the spacetime complex consists of just a collection of
  // vertices, presumably it can always be "oriented"?
  if (dimension < 1) return true;

  if (dimension == 1) {
    int i,vx[2];
    Graph* G = new Graph(nvertex);

    for(i=0; i<(signed) elements[1].size(); ++i) {
      elements[1][i].get_vertices(vx);
      G->add_edge(vx[0],vx[1]);
    }
    bool output = (G->bridge_count() == 0) ? true : false;
    delete G;
    return output;
  }

  if (nvertex == 1) return true;

  const int nf = (signed) elements[dimension-1].size(); 
  const int ns = (signed) elements[dimension].size();
  int i,j,m,k,o,in1,check,vx[1+dimension],orient[ns],sorient[nf],sorientn[nf];
  bool nzero,found,failed = false;
  hash_map face_index;
  std::vector<int>::const_iterator vit;
  std::stringstream s;
  std::vector<int> proc,nproc,vfacet,current;
  std::vector<int>* facets = new std::vector<int>[ns];

  for(i=0; i<ns; ++i) {
    elements[dimension][i].get_vertices(vx);
    for(j=0; j<=dimension; ++j) {
      vfacet.push_back(vx[j]);
    }
    facets[i] = vfacet;
    vfacet.clear();
  }

  for(i=0; i<nf; ++i) {
    elements[dimension-1][i].get_vertices(vx);
    for(j=0; j<dimension-1; ++j) {
      s << vx[j] << ":";
    }
    s << vx[dimension-1];
    face_index[s.str()] = i;
    s.str("");
  }
  
  for(i=0; i<ns; ++i) {
    orient[i] = 0;
  }
  orient[0] = 1;
  get_neighbours(ns,facets,0,proc);
  do {
    in1 = proc[0];
    current = facets[in1];
    nproc.clear();
    for(vit=proc.begin()+1; vit!=proc.end(); vit++) {
      nproc.push_back(*vit);
    }
    proc = nproc;
    if (orient[in1] != 0) continue;
    for(o=-1; o<2; o+=2) {
      induced_orientation(nf,current,o,face_index,sorientn);
      failed = false;
      for(i=0; i<ns; ++i) {
        if (orient[i] == 0) continue;
        m = 0;
        for(j=0; j<=dimension; ++j) {
          found = false;
          for(k=0; k<=dimension; ++k) {
            if (current[j] == facets[i][k]) {
              found = true;
              break;
            }
          }
          if (found) m++;
        }
        if (m < (dimension-1)) continue;
        induced_orientation(nf,facets[i],orient[i],face_index,sorient);
        for(j=0; j<nf; ++j) {
          check = sorientn[j] + sorient[j];
          if (check == 2 || check == -2) {
            failed = true;
            break;
          }
        }
        if (failed) break;
      }
      if (!failed) {
        orient[in1] = o;
        get_neighbours(ns,facets,in1,nproc);
        for(vit=nproc.begin(); vit!=nproc.end(); vit++) {
          proc.push_back(*vit);
        }
        break;
      }
    }
    nzero = true;
    for(i=0; i<nvertex; ++i) {
      if (orient[i] == 0) {
        nzero = false;
        break;
      }
    }
    if (nzero) break;
  } while(!failed && !proc.empty());
  delete[] facets;
  return !failed;
}

void Nexus::compute_homotopy(Group* pi1) const
{
  int i,j,q,r,ngen,n1,n2,n3,e1,e2,nf,vx[3],ntree;
  bool found;
  Word w(0);
  std::vector<int> tree_edges,generator,s1,s2;
  std::vector<Word> relations;
  const int ne = (signed) elements[1].size();
  const int nr = (dimension > 1) ? (signed) elements[2].size() : 0;

  pi1->clear();

  // First we need to calculate a spanning tree for the 1-skeleton of this complex...
  ntree = spanning_tree(tree_edges);
  for(i=0; i<ne; ++i) {
    s1.push_back(i);
  }
  for(i=0; i<nr; ++i) {
    s2.push_back(i);
  }

  // In princple, this is a group with ngenerators = nedges - ntree/2 and nrelations = ntriangles but we
  // can ignore those 2-simplices all of whose edges are in the spanning tree
  for(i=0; i<ne; ++i) {
    elements[1][s1[i]].get_vertices(vx);

    found = false;
    for(j=0; j<ntree; j+=2) {
      if (vx[0] == tree_edges[j] && vx[1] == tree_edges[j+1]) {
        found = true;
        break;
      }
    }

    if (!(found)) {
      // It's a generator...
      generator.push_back(vx[0]);
      generator.push_back(vx[1]);
    }
  }
  ngen = (signed) generator.size()/2;
  for(i=0; i<nr; ++i) {
    w.clear();
    elements[2][s2[i]].get_vertices(vx);

    nf = 0;

    e1 = vx[0];
    e2 = vx[1];
    found = false;
    n1 = -1;
    for(j=0; j<ngen; ++j) {
      if (e1 == generator[2*j] && e2 == generator[2*j+1]) {
        found = true;
        n1 = j;
        break;
      }
    }
    if (found) nf++;

    e1 = vx[0];
    e2 = vx[2];
    found = false;
    n2 = -1;
    for(j=0; j<ngen; ++j) {
      if (e1 == generator[2*j] && e2 == generator[2*j+1]) {
        found = true;
        n2 = j;
        break;
      }
    }
    if (found) nf++;

    e1 = vx[1];
    e2 = vx[2];
    found = false;
    n3 = -1;
    for(j=0; j<ngen; ++j) {
      if (e1 == generator[2*j] && e2 == generator[2*j+1]) {
        found = true;
        n3 = j;
        break;
      }
    }
    if (found) nf++;

    if (nf == 0) continue;

    if (nf == 1) {
      if (n1 >= 0) w.content.push_back(std::pair<int,int>(n1,1));
      if (n2 >= 0) w.content.push_back(std::pair<int,int>(n2,1));
      if (n3 >= 0) w.content.push_back(std::pair<int,int>(n3,1));
      relations.push_back(w);
      continue;
    }
    // At least two of the edges of this triangle are generators, we need to
    // order them correctly...
    if (n1 >= 0) {
      w.content.push_back(std::pair<int,int>(n1,1));
      q = generator[2*n1+1];

      if (n2 >= 0) {
        r = generator[2*n2];
        if (r == q) {
          w.content.push_back(std::pair<int,int>(n2,1));
          if (n3 >= 0) w.content.push_back(std::pair<int,int>(n3,1));
        }
        else {
          if (n3 >= 0) w.content.push_back(std::pair<int,int>(n3,1));
          w.content.push_back(std::pair<int,int>(n2,-1));
        }
        relations.push_back(w);
        continue;
      }
      r = generator[2*n3];
      if (r == q) {
        w.content.push_back(std::pair<int,int>(n3,1));
      }
      else {
        w.content.push_back(std::pair<int,int>(n3,-1));
      }
    }
    else {
      // So the two non-trivial generators are n2 and n3
      w.content.push_back(std::pair<int,int>(n2,1));
      q = generator[2*n2+1];
      r = generator[2*n3];
      if (q == r) {
        w.content.push_back(std::pair<int,int>(n3,1));
      }
      else {
        w.content.push_back(std::pair<int,int>(n3,-1));
      }
    }
    relations.push_back(w);
  }
  pi1->initialize(ngen,relations);
}
