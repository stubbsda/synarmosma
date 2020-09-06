#include "nexus.h"

using namespace SYNARMOSMA;

Nexus::Nexus() : Schema()
{

}

Nexus::Nexus(int D) : Schema()
{
  if (D < 1) throw std::invalid_argument("The dimensionality of a Nexus instance must be greater than zero!");

  elements = new std::vector<Cell>[1+D];
  index_table = new hash_map[1+D];
  topological_dimension = D;
}

Nexus::Nexus(const Nexus& source)
{
  clear();

  nvertex = source.nvertex;
  neighbours = source.neighbours;
  topological_dimension = source.topological_dimension;
  elements = new std::vector<Cell>[1+topological_dimension];
  index_table = new hash_map[1+topological_dimension];
  for(int i=0; i<=topological_dimension; ++i) {
    elements[i] = source.elements[i];
    index_table[i] = source.index_table[i];
  }
}

Nexus& Nexus::operator =(const Nexus& source)
{
  if (this == &source) return *this;

  clear();

  nvertex = source.nvertex;
  neighbours = source.neighbours;
  topological_dimension = source.topological_dimension;
  elements = new std::vector<Cell>[1+topological_dimension];
  index_table = new hash_map[1+topological_dimension];
  for(int i=0; i<=topological_dimension; ++i) {
    elements[i] = source.elements[i];
    index_table[i] = source.index_table[i];
  }

  return *this;
}

Nexus::~Nexus()
{
  if (topological_dimension > -1) {
    delete[] elements;
    delete[] index_table;
  }
}

void Nexus::initialize(int D)
{
  if (D < 1) throw std::invalid_argument("The dimensionality of a Nexus instance must be greater than zero!"); 

  clear();

  elements = new std::vector<Cell>[1+D];
  index_table = new hash_map[1+D];
  topological_dimension = D;
}

void Nexus::initialize(int D,int n)
{
  if (D < 1) throw std::invalid_argument("The dimensionality of a Nexus instance must be greater than zero!"); 
  if (n < 1) throw std::invalid_argument("The number of vertices in a Nexus instance must be greater than zero!"); 

  clear();

  elements = new std::vector<Cell>[1+D];
  index_table = new hash_map[1+D];
  topological_dimension = D;
  nvertex = n;
  std::set<int> null;
  for(int i=0; i<nvertex; ++i) {
    neighbours.push_back(null);
  }
}

void Nexus::clear()
{
  nvertex = 0;
  neighbours.clear();
  if (topological_dimension > -1) {
    delete[] elements;
    delete[] index_table;
  }
  topological_dimension = -1;
}

int Nexus::serialize(std::ofstream& s) const
{
  int i,j,n,count = 0;
  std::set<int>::const_iterator it;

  s.write((char*)(&nvertex),sizeof(int)); count += sizeof(int);
  s.write((char*)(&topological_dimension),sizeof(int)); count += sizeof(int);
  for(i=0; i<nvertex; ++i) {
    j = (signed) neighbours[i].size();
    s.write((char*)(&j),sizeof(int)); count += sizeof(int);
    for(it=neighbours[i].begin(); it!=neighbours[i].end(); ++it) {
      j = *it;
      s.write((char*)(&j),sizeof(int)); count += sizeof(int);
    }
  }
  for(i=1; i<topological_dimension; ++i) {
    n = (signed) elements[i].size();
    s.write((char*)(&n),sizeof(int)); count += sizeof(int);
    for(j=0; j<n; ++j) {
      count += elements[i][j].serialize(s);
    }
  }
  return count;
}

int Nexus::deserialize(std::ifstream& s)
{
  int i,j,k,n,count = 0;
  std::set<int> S;
  Cell sigma;

  clear();

  s.read((char*)(&nvertex),sizeof(int)); count += sizeof(int);
  s.read((char*)(&topological_dimension),sizeof(int)); count += sizeof(int);
  if (topological_dimension > -1) {
    elements = new std::vector<Cell>[1+topological_dimension];
    index_table = new hash_map[1+topological_dimension];
  }
  for(i=0; i<nvertex; ++i) {
    s.read((char*)(&n),sizeof(int)); count += sizeof(int);
    for(j=0; j<n; ++j) {
      s.read((char*)(&k),sizeof(int)); count += sizeof(int);
      S.insert(k);
    }
    neighbours.push_back(S);
    S.clear();
  }
  for(i=1; i<topological_dimension; ++i) {
    s.read((char*)(&n),sizeof(int)); count += sizeof(int);
    for(j=0; j<n; ++j) {
      count += sigma.deserialize(s);
      elements[i].push_back(sigma);
      index_table[i][sigma.vertices] = j;
    }
  }
  return count;
}

bool Nexus::paste(const Cell& c)
{
  std::set<int> vx; 
  int m = c.dimension();
  if (m > topological_dimension) throw std::invalid_argument("The Cell being pasted to this Nexus instance has too high a dimension!");

  c.get_vertices(vx);
  hash_map::const_iterator qt = index_table[m].find(vx);
  if (qt == index_table[m].end()) {
    elements[m].push_back(c);
    index_table[m][vx] = elements[m].size() - 1;
    // Now check if any vertices need to be added...
    int p = *vx.rbegin();
    if (p >= nvertex) {
      std::set<int> null;
      for(int i=nvertex; i<=p; ++i) {
        neighbours.push_back(null);
      }
      nvertex = 1 + p;
    }
    return true;
  }
  return false;
}

void Nexus::surface_construction(const std::string& type)
{
  // A method to construct certain standard surfaces for testing the correctness 
  // of the pseudomanifold and orientability routines...
  std::set<int> vx;
  
  std::string utype = boost::to_upper_copy(type);

  if (utype == "SPHERE") {
    // The 2-sphere $S^2$, which is an orientable pseudomanifold without boundary
    initialize(2,4);

    vx.insert(0); vx.insert(1); vx.insert(2);
    paste(vx);
    vx.clear();

    vx.insert(0); vx.insert(1); vx.insert(3);
    paste(vx);
    vx.clear();

    vx.insert(0); vx.insert(2); vx.insert(3);
    paste(vx);
    vx.clear();

    vx.insert(1); vx.insert(2); vx.insert(3);
    paste(vx);
    vx.clear();
  }
  else if (utype == "PROJECTIVE_PLANE") {
    // The real projective plane, a non-orientable pseudomanifold without boundary
    initialize(2,6);

    vx.insert(0); vx.insert(2); vx.insert(4);
    paste(vx);
    vx.clear();

    vx.insert(1); vx.insert(2); vx.insert(4);
    paste(vx);
    vx.clear();

    vx.insert(1); vx.insert(3); vx.insert(4);
    paste(vx);
    vx.clear();

    vx.insert(3); vx.insert(4); vx.insert(5);
    paste(vx);
    vx.clear();

    vx.insert(0); vx.insert(4); vx.insert(5);
    paste(vx);
    vx.clear();

    vx.insert(0); vx.insert(1); vx.insert(3);
    paste(vx);
    vx.clear();

    vx.insert(0); vx.insert(1); vx.insert(5);
    paste(vx);
    vx.clear();

    vx.insert(0); vx.insert(2); vx.insert(3);
    paste(vx);
    vx.clear();

    vx.insert(2); vx.insert(3); vx.insert(5);
    paste(vx);
    vx.clear();

    vx.insert(1); vx.insert(2); vx.insert(5);
    paste(vx);
    vx.clear();
  }
  else if (utype == "TORUS") {
    // The torus $S^1 \times S^1$, an orientable pseudomanifold without boundary
    initialize(2,9);

    vx.insert(0); vx.insert(1); vx.insert(3);
    paste(vx);
    vx.clear();

    vx.insert(1); vx.insert(3); vx.insert(4);
    paste(vx);
    vx.clear();

    vx.insert(1); vx.insert(2); vx.insert(4);
    paste(vx);
    vx.clear();

    vx.insert(2); vx.insert(4); vx.insert(5);
    paste(vx);
    vx.clear();

    vx.insert(0); vx.insert(2); vx.insert(5);
    paste(vx);
    vx.clear();

    vx.insert(0); vx.insert(3); vx.insert(5);
    paste(vx);
    vx.clear();

    vx.insert(3); vx.insert(4); vx.insert(6);
    paste(vx);
    vx.clear();

    vx.insert(4); vx.insert(6); vx.insert(7);
    paste(vx);
    vx.clear();

    vx.insert(4); vx.insert(5); vx.insert(7);
    paste(vx);
    vx.clear();

    vx.insert(5); vx.insert(7); vx.insert(8);
    paste(vx);
    vx.clear();

    vx.insert(3); vx.insert(5); vx.insert(8);
    paste(vx);
    vx.clear();

    vx.insert(3); vx.insert(6); vx.insert(8);
    paste(vx);
    vx.clear();

    vx.insert(0); vx.insert(6); vx.insert(7);
    paste(vx);
    vx.clear();

    vx.insert(0); vx.insert(1); vx.insert(7);
    paste(vx);
    vx.clear();

    vx.insert(1); vx.insert(7); vx.insert(8);
    paste(vx);
    vx.clear();

    vx.insert(1); vx.insert(2); vx.insert(8);
    paste(vx);
    vx.clear();

    vx.insert(2); vx.insert(6); vx.insert(8);
    paste(vx);
    vx.clear();

    vx.insert(0); vx.insert(2); vx.insert(6);
    paste(vx);
    vx.clear();
  }
  else if (utype == "MÖBIUS_STRIP") {
    // The Möbius strip, a non-orientable pseudomanifold with boundary
    initialize(2,8);

    vx.insert(0); vx.insert(1); vx.insert(2);
    paste(vx);
    vx.clear();

    vx.insert(1); vx.insert(2); vx.insert(3);
    paste(vx);
    vx.clear();

    vx.insert(2); vx.insert(3); vx.insert(4);
    paste(vx);
    vx.clear();

    vx.insert(3); vx.insert(4); vx.insert(5);
    paste(vx);
    vx.clear();

    vx.insert(4); vx.insert(5); vx.insert(6);
    paste(vx);
    vx.clear();

    vx.insert(5); vx.insert(6); vx.insert(7);
    paste(vx);
    vx.clear();

    vx.insert(6); vx.insert(7); vx.insert(1);
    paste(vx);
    vx.clear();

    vx.insert(0); vx.insert(1); vx.insert(7);
    paste(vx);
    vx.clear();
  }
  else {
    throw std::invalid_argument("Unrecognized surface type in Nexus::surface_construction method!");
  }
  regularization();
}

bool Nexus::consistent() const
{
  if (!Schema::consistent()) return false;
  // We shouldn't have any 0-simplices stored in these properties...
  if (!elements[0].empty() || !index_table[0].empty()) return false;

  int i,j,k,ns;
  std::set<int> vx;
  hash_map::const_iterator qt;

  for(i=topological_dimension; i>1; i--) {
    ns = (signed) elements[i].size();
    for(j=0; j<ns; ++j) {
      // Verify that every vertex in this simplex exists...
      elements[i][j].get_vertices(vx);
      k = *vx.rbegin();
      if (k >= nvertex) return false;
      // Now check the entailment property for an abstract simplicial complex...
      for(k=0; k<1+i; ++k) {
        qt = index_table[i-1].find(elements[i][j].faces[k]);
        if (qt == index_table[i-1].end()) return false;
      }
    }
  }

  return true;
}

int Nexus::compute_neighbour_graph(Graph* G) const
{
  if (topological_dimension < 1) throw std::runtime_error("Cannot compute the neighbour graph of a 0-dimensional simplicial complex!");

  int i,j;
  const int n = (signed) elements[topological_dimension].size();
  
  G->clear();

  // Loop over all d-simplices in the complex and add a vertex...
  for(i=0; i<n; ++i) {
    G->add_vertex();
  }

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i,j) schedule(dynamic,1)
#endif
  for(i=0; i<n; ++i) {
    for(j=1+i; j<n; ++j) {
      // Check if these two d-simplices are adjacent, i.e. if d of their respective 1+d vertices are the same.
      if (affinity(elements[topological_dimension][i],elements[topological_dimension][j]) == topological_dimension) {
#ifdef _OPENMP
#pragma omp critical
#endif
        G->add_edge(i,j);
      }
    }
  }
  return G->size();
}

void Nexus::compute_neighbours()
{
  if (nvertex != (signed) neighbours.size()) throw std::runtime_error("Illegal value of the nvertex property in Nexus::compute_neighbours!");
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
  hash_map::const_iterator qt;

  for(i=topological_dimension; i>1; i--) {
    ns = (signed) elements[i].size();
    for(j=0; j<ns; ++j) {
      for(k=0; k<1+i; ++k) {
        qt = index_table[i-1].find(elements[i][j].faces[k]);
        if (qt == index_table[i-1].end()) throw std::runtime_error("Missing entourage element in Nexus::compute_entourages!");
        elements[i-1][qt->second].entourage.insert(j);
      }
    }
  }
}

int Nexus::regularization()
{
  // If a given d-simplex exists, we need to make sure that all of
  // its sub-simplices also exist in the nexus
  int i,j,k,n,m,na = 0;
  Cell S;
  hash_map::const_iterator qt;

  for(i=topological_dimension; i>=2; i--) {
    n = (signed) elements[i].size();
    m = (signed) elements[i-1].size();
    for(j=0; j<n; ++j) {
      for(k=0; k<1+i; ++k) {
        qt = index_table[i-1].find(elements[i][j].faces[k]);
        if (qt == index_table[i-1].end()) {
          S = Cell(elements[i][j].faces[k]);
          elements[i-1].push_back(S);
          index_table[i-1][S.vertices] = m;
          m++;
          na++;
        }
      }
    }
  }
  compute_entourages();
  compute_neighbours();

  if (!consistent()) throw std::runtime_error("The Nexus instance is inconsistent after regularization!");

  return na;
}

void Nexus::closure(const std::set<std::set<int> >& S,Nexus* NX,int* offset) const
{
  int i,n,m;
  unsigned int j;
  std::set<int> vx,s;
  std::vector<Cell>* cell_array = new std::vector<Cell>[1+topological_dimension];
  hash_map::const_iterator qt;
  std::set<std::set<int> >::const_iterator st;
  std::set<int>::const_iterator it;

  for(st=S.begin(); st!=S.end(); st++) {
    s = *st;
    if (s.empty()) continue;
    n = s.size() - 1;
    if (n > 0) {
      qt = index_table[n].find(s);
      m = qt->second;
      cell_array[n].push_back(elements[n][m]);
    }
    else {
      vx.insert(*(s.begin()));
    }
  }

  for(i=0; i<nvertex; ++i) {
    offset[i] = -1;
  }
  // Here n will be the number of vertices in the closure of the new Nexus while m 
  // will be its dimension
  n = 0;
  m = 0;
  for(i=1; i<=topological_dimension; ++i) {
    if (cell_array[i].empty()) continue;
    m = i;
    for(j=0; j<cell_array[i].size(); ++j) {
      for(it=cell_array[i][j].vertices.begin(); it!=cell_array[i][j].vertices.end(); ++it) {
        if (offset[*it] == -1) {
          offset[*it] = n;
          n++;
        }
      }
    }
  }
  for(it=vx.begin(); it!=vx.end(); ++it) {
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
      for(it=cell_array[i][j].vertices.begin(); it!=cell_array[i][j].vertices.end(); ++it) {
        vx.insert(offset[*it]);
      }
      NX->paste(vx);
      vx.clear();
    }
  }
  NX->regularization();
  delete[] cell_array;
}

void Nexus::link(const std::set<std::set<int> >& S,std::vector<Cell>* cell_array) const
{
  int i,offset1[nvertex],offset2[nvertex];
  unsigned int j,k;
  bool found;
  std::vector<Cell>* cell_array1 = new std::vector<Cell>[1+topological_dimension];
  std::vector<Cell>* cell_array2 = new std::vector<Cell>[1+topological_dimension];
  std::set<int> vx;
  std::set<std::set<int> > sset;
  std::set<int>::const_iterator it;
  Nexus* N1 = new Nexus;
  Nexus* N2 = new Nexus;

  // The link of S is defined to be closure(star(S)) - star(closure(S)), we begin 
  // with the first term
  star(S,cell_array1);
  // Now we need to convert cell_array to a set of strings...
  for(i=0; i<=topological_dimension; ++i) {
    for(j=0; j<cell_array1[i].size(); ++j) {
      sset.insert(cell_array1[i][j].vertices);
    }
  }  
  closure(sset,N1,offset1);

  // Now the second term...
  closure(S,N2,offset2);
  sset.clear();
  for(i=0; i<=N2->topological_dimension; ++i) {
    for(j=0; j<N2->elements[i].size(); ++j) {
      for(it=N2->elements[i][j].vertices.begin(); it!=N2->elements[i][j].vertices.end(); ++it) {
        vx.insert(offset2[*it]);
      }
      sset.insert(vx);
      vx.clear();
    }
  }
  star(sset,cell_array2);

  // Finally, construct the output by taking every cell that exists in N1 and, if it doesn't also 
  // exist in cell_array2, adding it to cell_array
  for(i=0; i<=N1->topological_dimension; ++i) {
    for(j=0; j<N1->elements[i].size(); ++j) {
      for(it=N1->elements[i][j].vertices.begin(); it!=N1->elements[i][j].vertices.end(); ++it) {
        vx.insert(offset1[*it]);
      }
      // Does this element exist somewhere in cell_array2?
      found = false;
      for(k=0; k<cell_array2[i].size(); ++k) {
        if (vx == cell_array2[i][k].vertices) {
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

void Nexus::star(const std::set<std::set<int> >& S,std::vector<Cell>* output) const
{
  int n,m,d,nc;
  std::set<int> vx;
  hash_map ckey;
  std::vector<Cell> clist;
  hash_map::const_iterator qt;
  std::set<std::set<int> >::const_iterator it;

  for(it=S.begin(); it!=S.end(); ++it) {
    vx = *it;
    n = vx.size() - 1; 
    if (n > 0) {
      qt = index_table[n-1].find(vx);
      m = qt->second;
    }
    else {
      m = *(vx.begin());
    }
    ascend(n,m,clist);
  }
  for(n=0; n<=topological_dimension; ++n) {
    output[n].clear();
  }
  nc = (signed) clist.size();
  d = clist[0].dimension();
  output[d].push_back(clist[0]);
  ckey[clist[0].vertices] = 0;
  m = 1;
  for(n=1; n<nc; ++n) {
    qt = ckey.find(clist[n].vertices);
    if (qt == ckey.end()) {
      d = clist[n].dimension();
      output[d].push_back(clist[n]);
      ckey[clist[n].vertices] = m;
      m++;
    }
  }
}

void Nexus::ascend(int D,int n,std::vector<Cell>& out) const
{
  if (D < 1 && D > topological_dimension) return;

  int i;
  std::set<int>::const_iterator it;
  std::set<int> S = elements[D][n].entourage;

  for(it=S.begin(); it!=S.end(); ++it) {
    i = *it;
    out.push_back(elements[D+1][i]);
    ascend(1+D,i,out);
  }
}

bool Nexus::pseudomanifold(bool* boundary) const
{
  int i,j;
  bool found;
  std::vector<int> candidates;

  if (topological_dimension < 1) return false;

  *boundary = false;

  // First the non-branching requirement...
  for(i=0; i<(signed) elements[topological_dimension-1].size(); ++i) {
    j = (signed) elements[topological_dimension-1][i].entourage.size();
    // If the entourage size is ever 1, this is a pseudomanifold with boundary!
    if (j > 2 || j < 1) return false;
    if (j == 1) *boundary = true;
  }

  // Next, check the dimensional homogeneity...
  const int nd = (signed) elements[topological_dimension].size();

  for(i=0; i<nvertex; ++i) {
    found = false;
    for(j=0; j<nd; ++j) {
      if (elements[topological_dimension][j].contains(i)) {
        found = true;
        break;
      }
    }
    if (!found) return false;
  }

  // Finally, as we have verified that this is a non-branching simplicial complex, 
  // we next need to check if it is strongly connected. 
  // If there's only a single d-simplex, then of course it's strong connected...
  if (nd == 1) return true;
  Graph G;
  compute_neighbour_graph(&G);
  // If the graph of d-simplices based on their adjacency is connected, then the 
  // complex is strongly connected. 
  return G.connected();
}

bool Nexus::orientable() const
{
  // Program to calculate whether the complex is orientable, note
  // this only makes sense if the complex is also a pseudomanifold!

  // If the spacetime complex consists of just a collection of
  // vertices, presumably it can always be "oriented"?
  if (topological_dimension < 1) return true;

  if (topological_dimension == 1) {
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

  const int nf = (signed) elements[topological_dimension-1].size(); 
  const int ns = (signed) elements[topological_dimension].size();

  if (ns == 1) return true;

  int i,j,m,k,o,in1,check,vx[1+topological_dimension],orient[ns],sorient[nf],sorientn[nf];
  bool nzero,found,failed = false;
  hash_map face_index;
  std::vector<int>::const_iterator vit;
  std::set<int> S;
  std::vector<int> proc,nproc,vfacet,current;
  std::vector<int>* facets = new std::vector<int>[ns];

  for(i=0; i<ns; ++i) {
    elements[topological_dimension][i].get_vertices(vx);
    for(j=0; j<=topological_dimension; ++j) {
      vfacet.push_back(vx[j]);
    }
    facets[i] = vfacet;
    vfacet.clear();
  }

  for(i=0; i<nf; ++i) {
    elements[topological_dimension-1][i].get_vertices(S);
    face_index[S] = i;
  }
  
  for(i=0; i<ns; ++i) {
    orient[i] = 0;
  }
  orient[0] = 1;
  assemble_neighbours(ns,facets,0,proc);
  if (proc.empty()) throw std::runtime_error("Error in neighbour assembly for Nexus::orientable method!");
  do {
    in1 = proc[0];
    current = facets[in1];
    nproc.clear();
    for(vit=proc.begin()+1; vit!=proc.end(); ++vit) {
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
        for(j=0; j<=topological_dimension; ++j) {
          found = false;
          for(k=0; k<=topological_dimension; ++k) {
            if (current[j] == facets[i][k]) {
              found = true;
              break;
            }
          }
          if (found) m++;
        }
        if (m < (topological_dimension-1)) continue;
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
        assemble_neighbours(ns,facets,in1,nproc);
        for(vit=nproc.begin(); vit!=nproc.end(); ++vit) {
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

