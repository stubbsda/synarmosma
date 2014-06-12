#include "geometry.h"

extern Random RND;

const int Geometry::background_dimension;

Geometry::Geometry()
{
  nvertex = 0;
  euclidean = true;
  relational = false;
  uniform = true;
  vperturb = -1;
}

Geometry::Geometry(bool type,bool model,bool flat)
{
  nvertex = 0;
  euclidean = type;
  relational = model;
  uniform = flat;
  vperturb = -1;
}

Geometry::Geometry(const Geometry& source)
{
  nvertex = source.nvertex;
  euclidean = source.euclidean;
  relational = source.relational;
  uniform = source.uniform;
  vperturb = source.vperturb;
  original = source.original;
  distances = source.distances;
  index = source.index;
  coordinates = source.coordinates;
}

Geometry& Geometry::operator =(const Geometry& source)
{
  if (this == &source) return *this;

  nvertex = source.nvertex;
  euclidean = source.euclidean;
  relational = source.relational;
  uniform = source.uniform;
  vperturb = source.vperturb;
  original = source.original;
  distances = source.distances;
  index = source.index;
  coordinates = source.coordinates;

  return *this;
}

Geometry::~Geometry()
{
  clear();
}

void Geometry::clear()
{
  if (relational) {
    index.clear();
  }
  else {
    for(int i=0; i<nvertex; ++i) {
      coordinates[i].clear();
    }
    coordinates.clear();
  }
  distances.clear();
  original.clear();
  vperturb = -1;
  nvertex = 0;
  euclidean = true;
  relational = false;
  uniform = true;
}

void Geometry::serialize(std::ofstream& s) const
{
  int i,j;
  double x;

  s.write((char*)(&nvertex),sizeof(int));
  s.write((char*)(&euclidean),sizeof(bool));
  s.write((char*)(&relational),sizeof(bool));
  s.write((char*)(&uniform),sizeof(bool));
  if (relational) {
    hash_map::const_iterator qt;
    for(i=0; i<nvertex; ++i) {
      for(j=1+i; j<nvertex; ++j) {
        qt = index.find(make_key(i,j));
        x = distances[qt->second];
        s.write((char*)(&x),sizeof(double));
      }
    }
  }
  else {
    int n;
    for(i=0; i<nvertex; ++i) {
      n = (signed) coordinates[i].size();
      s.write((char*)(&n),sizeof(int));
      for(j=0; j<n; ++j) {
        x = coordinates[i][j];
        s.write((char*)(&x),sizeof(double));
      }
    }
  }
}

void Geometry::deserialize(std::ifstream& s)
{
  int i,j,n;
  double x;

  clear();

  s.read((char*)(&nvertex),sizeof(int));
  s.read((char*)(&euclidean),sizeof(bool));
  s.read((char*)(&relational),sizeof(bool));
  s.read((char*)(&uniform),sizeof(bool));
  if (relational) {
    n = 0;
    for(i=0; i<nvertex; ++i) {
      for(j=1+i; j<nvertex; ++j) {
        s.read((char*)(&x),sizeof(double));
        index[make_key(i,j)] = n;
        distances.push_back(x);
        n++;
      }
    }
  }
  else {
    std::vector<double> xc;
    for(i=0; i<nvertex; ++i) {
      s.read((char*)(&n),sizeof(int));
      for(j=0; j<n; ++j) {
        s.read((char*)(&x),sizeof(double));
        xc.push_back(x);
      }
      coordinates.push_back(xc);
      xc.clear();
    }
  }
  compute_distances();
}

double Geometry::dot_product(const std::vector<double>& vx,const std::vector<double>& vy) const
{
  double output = 0.0;
  if (uniform) {
    for(int i=0; i<Geometry::background_dimension; ++i) {
      output += vx[i]*vy[i];
    }
    return output;
  }
  int i,k = 0,n = (signed) vx.size(),m = (signed) vy.size();
  std::vector<double> vlx = vx;
  std::vector<double> vly = vy;
  int l = n;

  if (n < m) {
    for(i=n; i<m; ++i) {
      vlx.push_back(RND.nrandom(0.0,1.0/double(2+k)));
      k++;
    }
    l = m;
  }
  else if (m < n) {
    for(i=m; i<n; ++i) {
      vly.push_back(RND.nrandom(0.0,1.0/double(2+k)));
      k++;
    }
  }
  for(i=0; i<l; ++i) {
    output += vlx[i]*vly[i];
  }
  return output;
}

void Geometry::get_implied_vertices(int n,std::set<int>& vx) const
{
  int i,j;
  if (relational) {
    // Loop over all vertex pairs and find the one that corresponds to index
    // "n"
    hash_map::const_iterator qt;
    for(i=0; i<nvertex; ++i) {
      for(j=1+i; j<nvertex; ++j) {
        qt = index.find(make_key(i,j));
        if (qt->second == n) {
          vx.insert(i);
          vx.insert(j);
          return;
        }
      }
    }
  }
  else {
    int kt = 0;
    for(i=0; i<nvertex; ++i) {
      for(j=0; j<(signed) coordinates[i].size(); ++j) {
        if (kt == n) {
          vx.insert(i);
          return;
        }
        kt++;
      }
    }
  }
}

void Geometry::load(const Geometry* source)
{
  nvertex = source->nvertex;
  euclidean = source->euclidean;
  distances = source->distances;
  if (relational) {
    index = source->index;
  }
  else {
    coordinates = source->coordinates;
  }
}

void Geometry::store(Geometry* target) const
{
  target->nvertex = nvertex;
  target->euclidean = euclidean;
  target->distances = distances;
  if (relational) {
    target->index = index;
  }
  else {
    target->coordinates = coordinates;
  }
}

int Geometry::vertex_order(int n,int m) const
{
  if (relational || euclidean) return -1;

  double sum = -(coordinates[n][0] - coordinates[m][0])*(coordinates[n][0] - coordinates[m][0]);
  for(int i=1; i<Geometry::background_dimension; ++i) {
    sum += (coordinates[n][i] - coordinates[m][i])*(coordinates[n][i] - coordinates[m][i]);
  }
  if (sum > 0.0) return -1;
  int output = (coordinates[n][0] < coordinates[m][0]) ? 1 : 0;
  return output;
}

void Geometry::initialize(int n,const std::string& type)
{
  if (!relational) {
    std::cerr << "Illegal geometric method call (initialize) for absolute model!" << std::endl;
    std::exit(1);
  }

  distances.clear();
  index.clear();

  int i,j,k = 0;
  if (type == "CARTESIAN") {
    int l,p,q,s,alpha;
    double delta;
    std::vector<int> x,y;

    nvertex = ipow(n,Geometry::background_dimension);
    for(i=0; i<nvertex; ++i) {
      s = i;
      for(j=Geometry::background_dimension; j>=1; j--) {
        p = ipow(n,j-1);
        q = s/p;
        x.push_back(q);
        s -= p*q;
      }
      x.push_back(s);
      for(j=1+i; j<nvertex; ++j) {
        s = j;
        for(l=Geometry::background_dimension; l>=1; --l) {
          p = ipow(n,l-1);
          q = s/p;
          y.push_back(q);
          s -= p*q;
        }
        y.push_back(s);
        delta = double((x[0] - y[0])*(x[0] - y[0]));
        if (!euclidean) delta = -delta;
        for(l=1; l<Geometry::background_dimension; ++l) {
          alpha = (x[l] - y[l])*(x[l] - y[l]);
          delta += double(alpha);
        }
        distances.push_back(delta);
        index[make_key(i,j)] = k;
        k++;
        y.clear();
      }
      x.clear();
    }
  }
  else if (type == "SINGLETON") {
    // Nothing to do here as there's no relational geometry with just a
    // single vertex...
    nvertex = 1;
  }
  else if (type == "RANDOM") {
    nvertex = n;
    if (euclidean) {
      for(i=0; i<nvertex; ++i) {
        for(j=1+i; j<nvertex; ++j) {
          distances.push_back(0.25 + 10.0*RND.drandom());
          index[make_key(i,j)] = k;
          k++;
        }
      }
    }
    else {
      double alpha,r = 1.0 + double(Geometry::background_dimension-1);
      for(i=0; i<nvertex; ++i) {
        for(j=1+i; j<nvertex; ++j) {
          alpha = -5.0 + 5.0*r*RND.drandom();
          distances.push_back(alpha);
          index[make_key(i,j)] = k;
          k++;
        }
      }
    }
  }
  else if (type == "MONOPLEX") {
    // In this case n is the monoplex dimension, so the number of vertices
    // is 1+n
    nvertex = 1 + n;
    if (euclidean) {
      for(i=1; i<nvertex; ++i) {
        distances.push_back(1.0);
        index[make_key(0,i)] = k;
        k++;
        for(j=1+i; j<nvertex; ++j) {
          distances.push_back(2.0);
          index[make_key(i,j)] = k;
          k++;
        }
      }
    }
    else {
      distances.push_back(-1.0);
      index["0:1"] = k;
      k++;
      for(i=2; i<nvertex; ++i) {
        distances.push_back(1.0);
        index[make_key(0,i)] = k;
        k++;
        distances.push_back(0.0);
        index[make_key(1,i)] = k;
        k++;
        for(j=1+i; j<nvertex; ++j) {
          distances.push_back(2.0);
          index[make_key(i,j)] = k;
          k++;
        }
      }
    }
  }
  else {
    std::cerr << "Illegal spacetime type in Geometry::initialize!" << std::endl;
    std::exit(1);
  }
}

void Geometry::vertex_difference(int n,int m,std::vector<double>& delta) const
{
  if (relational) {
    std::cerr << "Illegal geometric method call (vertex_difference) for relational model!" << std::endl;
    std::exit(1);
  }
  int i;

  delta.clear();
  if (uniform) {
    for(i=0; i<Geometry::background_dimension; ++i) {
      delta.push_back(coordinates[n][i] - coordinates[m][i]);
    }
    return;
  }
  int d1 = (signed) coordinates[n].size();
  int d2 = (signed) coordinates[m].size();
  int k=0,l = d1;

  std::vector<double> vlx = coordinates[n];
  std::vector<double> vly = coordinates[m];

  if (d1 < d2) {
    for(i=d1; i<d2; ++i) {
      vlx.push_back(RND.nrandom(0.0,1.0/double(2+k)));
      ++k;
    }
    l = d2;
  }
  else if (d2 < d1) {
    for(i=d2; i<d1; ++i) {
      vly.push_back(RND.nrandom(0.0,1.0/double(2+k)));
      ++k;
    }
  }

  for(i=0; i<l; ++i) {
    delta.push_back(vlx[i] - vly[i]);
  }
}

double Geometry::inner_product(const std::vector<double>* laplacian,const std::vector<int>& offset,int nv) const
{
  if (relational) return 0.0;

  int i,j,k,l,nelements,column;
  double value,result[nv][Geometry::background_dimension],sum[Geometry::background_dimension];
  double energy = 0.0;

  for(i=0; i<nv; ++i) {
    for(j=0; j<Geometry::background_dimension; ++j) {
      sum[j] = 0.0;
    }
    l = 0;
    nelements = (signed) laplacian[i].size()/2;
    for(j=0; j<nelements; ++j) {
      value = laplacian[i][l];
      column = int(laplacian[i][l+1]);
      for(k=0; k<Geometry::background_dimension; ++k) {
        sum[k] += coordinates[column][k]*value;
      }
      l += 2;
    }
    for(j=0; j<Geometry::background_dimension; ++j) {
      result[i][j] = sum[j];
    }
  }

  for(i=0; i<Geometry::background_dimension; ++i) {
    value = 0.0;
    for(j=0; j<nvertex; ++j) {
      if (offset[j] == -1) continue;
      value += result[offset[j]][i]*coordinates[j][i];
    }
    energy += value;
  }
  return energy;
}

void Geometry::reciprocate()
{
  if (relational) {
    const int n = (signed) distances.size();
    for(int i=0; i<n; ++i) {
      distances[i] = 1.0/distances[i];
    }
  }
  else {

  }
}

void Geometry::rollback()
{
  if (relational) {
    int i,j = 0;
    hash_map::const_iterator qt;
    for(i=0; i<nvertex; ++i) {
      if (i == vperturb) continue;
      qt = index.find(make_key(i,vperturb));
      distances[qt->second] = original[j];
      j++;
    }
  }
  else {
    coordinates[vperturb] = original;
  }
}

void Geometry::perturb_vertex(int v)
{
  vperturb = v;
  if (relational) {
    int i,j;
    hash_map::const_iterator qt;

    original.clear();
    for(i=0; i<nvertex; ++i) {
      if (i == v) continue;
      qt = index.find(make_key(i,v));
      j = qt->second;
      original.push_back(distances[j]);
      distances[j] += RND.nrandom(0.0,0.5);
    }
  }
  else {
    int k = RND.irandom(Geometry::background_dimension);
    original = coordinates[v];
    coordinates[v][k] += RND.nrandom(0.0,0.1);
  }
}

void Geometry::add_vertex(const std::set<int>& antecedents)
{
  if (antecedents.empty()) {
    add_vertex(-1);
  }
  else {
    if (RND.drandom() < 0.25) {
      add_vertex(RND.irandom(antecedents));
    }
    else {
      int i,j;
      double l,na = double(antecedents.size());
      std::set<int>::const_iterator it;
      if (relational) {
        hash_map::const_iterator qt;

        j = (signed) distances.size();
        for(i=0; i<nvertex; ++i) {
          l = 0.0;
          for(it=antecedents.begin(); it!=antecedents.end(); it++) {
            qt = index.find(make_key(i,*it));
            l += distances[qt->second];
          }
          l = l/na;
          index[make_key(i,nvertex)] = j;
          distances.push_back(l);
          j++;
        }
        nvertex++;
      }
      else {
        std::vector<double> xc,avg_x;
        for(i=0; i<Geometry::background_dimension; ++i) {
          avg_x.push_back(0.0);
        }
        for(it=antecedents.begin(); it!=antecedents.end(); it++) {
          j = *it;
          for(i=0; i<Geometry::background_dimension; ++i) {
            avg_x[i] += coordinates[j][i];
          }
        }
        for(i=0; i<Geometry::background_dimension; ++i) {
          avg_x[i] /= na;
        }
        for(i=0; i<Geometry::background_dimension; ++i) {
          xc.push_back(RND.nrandom(avg_x[i],0.5));
        }
        add_vertex(xc);
      }
    }
  }
}

void Geometry::add_vertex(int parent,double mutation)
{
  int i,j;
  double alpha;
  std::vector<double> x;

  if (parent == -1) {
    // No antecedent, so place the vertex randomly...
    if (relational) {
      j = (signed) distances.size();
      if (euclidean) {
        for(i=0; i<nvertex; ++i) {
          alpha = 1.0 + 15.0*RND.drandom();
          distances.push_back(alpha);
          index[make_key(i,nvertex)] = j;
          j++;
        }
      }
      else {
        // Skew the random number generation so the inter-vertex
        // lengths are properly weighted between timelike and
        // spacelike...
        double r = 1.0 + double(Geometry::background_dimension-1);
        for(i=0; i<nvertex; ++i) {
          alpha = -5.0 + 5.0*r*RND.drandom();
          distances.push_back(alpha);
          index[make_key(i,nvertex)] = j;
          j++;
        }
      }
    }
    else {
      for(i=0; i<Geometry::background_dimension; ++i) {
        alpha = -10.0 + 20.0*RND.drandom();
        x.push_back(alpha);
      }
      coordinates.push_back(x);
    }
  }
  else {
    // Should be close to its parent vertex...
    if (relational) {
      hash_map::const_iterator qt;
      j = (signed) distances.size();
      for(i=0; i<nvertex; ++i) {
        if (i == parent) continue;
        qt = index.find(make_key(i,parent));
        alpha = RND.nrandom(distances[qt->second],mutation/10.0);
        distances.push_back(alpha);
        index[make_key(i,nvertex)] = j;
        j++;
      }
      alpha = RND.nrandom(0.0,0.1);
      distances.push_back(alpha);
      index[make_key(parent,nvertex)] = j;
    }
    else {
      int q,p = RND.irandom(Geometry::background_dimension);
      double r = RND.drandom(0.1+0.5*mutation,0.2+mutation);
      x = coordinates[parent];
      alpha = RND.drandom(0.0,2.0*M_PI);
      do {
        q = RND.irandom(Geometry::background_dimension);
        if (q != p) break;
      } while(true);
      x[p] += r*std::cos(alpha);
      x[q] += r*std::sin(alpha);
      coordinates.push_back(x);
    }
  }
  nvertex++;
}

void Geometry::geometry_modification(int n,double mu,double sigma)
{
  vperturb = n;
  original.clear();
  if (relational) {
    original.push_back(distances[n]);
    distances[n] *= RND.nrandom(mu,sigma);
    return;
  }
  int i,j,kt = 0;
  for(i=0; i<nvertex; ++i) {
    for(j=0; j<(signed) coordinates[i].size(); ++j) {
      if (kt == n) {
        original.push_back(coordinates[i][j]);
        coordinates[i][j] *= RND.nrandom(mu,sigma);
        return;
      }
      kt++;
    }
  }
}

void Geometry::geometry_restoration()
{
  if (relational) {
    distances[vperturb] = original[0];
    return;
  }
  int i,j,kt = 0;
  for(i=0; i<nvertex; ++i) {
    for(j=0; j<(signed) coordinates[i].size(); ++j) {
      if (kt == vperturb) {
        coordinates[i][j] = original[0];
        return;
      }
      kt++;
    }
  }
}

void Geometry::multiplicative_modification(int v,bool total,double mu,double sigma)
{
  int i;

  vperturb = v;
  if (relational) {
    int j;
    hash_map::const_iterator qt;

    original.clear();
    for(i=0; i<nvertex; ++i) {
      if (i == v) continue;
      qt = index.find(make_key(v,i));
      j = qt->second;
      original.push_back(distances[j]);
      distances[j] *= RND.nrandom(mu,sigma);
    }
  }
  else {
    original = coordinates[v];
    if (total) {
      for(i=0; i<(signed) coordinates[v].size(); ++i) {
        coordinates[v][i] *= RND.nrandom(mu,sigma);
      }
    }
    else {
      i = (signed) coordinates[v].size();
      coordinates[v][RND.irandom(i)] *= RND.nrandom(mu,sigma);
    }
  }
}

void Geometry::additive_modification(int v,bool total,double mu,double sigma)
{
  int i;

  vperturb = v;
  if (relational) {
    int j;
    hash_map::const_iterator qt;

    original.clear();
    for(i=0; i<nvertex; ++i) {
      if (i == v) continue;
      qt = index.find(make_key(v,i));
      j = qt->second;
      original.push_back(distances[j]);
      distances[j] += RND.nrandom(mu,sigma);
    }
  }
  else {
    original = coordinates[v];
    if (total) {
      for(i=0; i<(signed) coordinates[v].size(); ++i) {
        coordinates[v][i] += RND.nrandom(mu,sigma);
      }
    }
    else {
      i = (signed) coordinates[v].size();
      coordinates[v][RND.irandom(i)] += RND.nrandom(mu,sigma);
    }
  }
}

bool Geometry::adjust_dimension(const std::vector<int>& vdimension)
{
  if (relational) return false;
  
  int i,j,n,m;
  std::set<int> vmodified;
  std::vector<double> x;

  for(i=0; i<nvertex; ++i) {
    if (vdimension[i] == -1) continue;
    // First check to see if the dimension has changed...
    m = vdimension[i];
    n = (signed) coordinates[i].size();
    if (m <= Geometry::background_dimension) {
      if (n != Geometry::background_dimension) {
        vmodified.insert(i);
        x = coordinates[i];
        coordinates[i].clear();
        for(j=0; j<Geometry::background_dimension; ++j) {
          coordinates[i].push_back(x[j]);
        }
      }
    }
    else {
      if (!uniform) {
        if (n != m) {
          vmodified.insert(i);
          if (n > m) {
            x = coordinates[i];
            coordinates[i].clear();
            for(j=0; j<m; ++j) {
              coordinates[i].push_back(x[j]);
            }
          }
          else {
            for(j=n; j<m; ++j) {
              coordinates[i].push_back(RND.nrandom(0.0,1.0));
            }
          }
        }
      }
    }
  }
  if (!vmodified.empty()) {
    compute_distances(vmodified);
    return true;
  }
  return false;
}

void Geometry::compute_distances(const std::set<int>& vmodified)
{
  if (relational) return;

  int i,j,k,in1;
  double delta;
  std::set<int>::const_iterator it;

  if (euclidean) {
    for(it=vmodified.begin(); it!=vmodified.end(); it++) {
      i = *it;
      for(j=0; j<i; ++j) {
        in1 = j*nvertex - j*(j+1)/2;
        delta = 0.0;
        for(k=0; k<Geometry::background_dimension; ++k) {
          delta += (coordinates[i][k] - coordinates[j][k])*(coordinates[i][k] - coordinates[j][k]);
        }
        distances[in1+i-(1+j)] = delta;
      }
      in1 = i*nvertex - i*(i+1)/2;
      for(j=1+i; j<nvertex; ++j) {
        delta = 0.0;
        for(k=0; k<Geometry::background_dimension; ++k) {
          delta += (coordinates[i][k] - coordinates[j][k])*(coordinates[i][k] - coordinates[j][k]);
        }
        distances[in1+j-(1+i)] = delta;
      }
    }
  }
  else {
    for(it=vmodified.begin(); it!=vmodified.end(); it++) {
      i = *it;
      for(j=0; j<i; ++j) {
        in1 = j*nvertex - j*(j+1)/2;
        delta = -(coordinates[i][0] - coordinates[j][0])*(coordinates[i][0] - coordinates[j][0]);
        for(k=1; k<Geometry::background_dimension; ++k) {
          delta += (coordinates[i][k] - coordinates[j][k])*(coordinates[i][k] - coordinates[j][k]);
        }
        distances[in1+i-(1+j)] = delta;
      }
      in1 = i*nvertex - i*(i+1)/2;
      for(j=1+i; j<nvertex; ++j) {
        delta = -(coordinates[i][0] - coordinates[j][0])*(coordinates[i][0] - coordinates[j][0]);
        for(k=1; k<Geometry::background_dimension; ++k) {
          delta += (coordinates[i][k] - coordinates[j][k])*(coordinates[i][k] - coordinates[j][k]);
        }
        distances[in1+j-(1+i)] = delta;
      }
    }
  }
}

void Geometry::compute_distances()
{
  if (relational) return;

  int i,j,k,in1;
  double delta;

  distances.clear();
  for(i=0; i<nvertex*(nvertex-1)/2; ++i) {
    distances.push_back(0.0);
  }
  if (euclidean) {
#ifdef PARALLEL
    #pragma omp parallel for default(shared) private(i,j,k,in1,delta) schedule(dynamic,1)
#endif
    for(i=0; i<nvertex; ++i) {
      in1 = i*nvertex - i*(i+1)/2;
      for(j=1+i; j<nvertex; ++j) {
        delta = 0.0;
        for(k=0; k<Geometry::background_dimension; ++k) {
          delta += (coordinates[i][k] - coordinates[j][k])*(coordinates[i][k] - coordinates[j][k]);
        }
        distances[in1+j-(1+i)] = delta;
      }
    }
  }
  else {
#ifdef PARALLEL
    #pragma omp parallel for default(shared) private(i,j,k,in1,delta) schedule(dynamic,1)
#endif
    for(i=0; i<nvertex; ++i) {
      in1 = i*nvertex - i*(i+1)/2;
      for(j=1+i; j<nvertex; ++j) {
        delta = -(coordinates[i][0] - coordinates[j][0])*(coordinates[i][0] - coordinates[j][0]);
        for(k=1; k<Geometry::background_dimension; ++k) {
          delta += (coordinates[i][k] - coordinates[j][k])*(coordinates[i][k] - coordinates[j][k]);
        }
        distances[in1+j-(1+i)] = delta;
      }
    }
  }
}

int Geometry::compute_coordinates(std::vector<double>& x) const
{
  int i,j,k,edim = 0;
  x.clear();

  if (!relational) {
    // Calculate the maximum simplicial dimension...
    for(i=0; i<nvertex; ++i) {
      j = (signed) coordinates[i].size();
      if (j > edim) edim = j;
    }
    // Pad with zeroes as necessary...
    for(i=0; i<nvertex; ++i) {
      k = (signed) coordinates[i].size();
      for(j=0; j<k; ++j) {
        x.push_back(coordinates[i][j]);
      }
      for(j=k; j<edim; ++j) {
        x.push_back(0.0);
      }
    }
    return edim;
  }

  assert(euclidean);
  /*
  int k,in1,its = 0;
  double b,lambda,err,err_old = 0.0,delta,sum;
  hash_map::const_iterator qt;
  std::vector<double> xnew,cdistance,BZ,Delta;
  const int M = 150;
  const double cutoff = 0.00001;
  const double pfactor = 1.0/double(nvertex);

  for(i=0; i<nvertex; ++i) {
    in1 = nvertex*i - i*(i+1)/2;
    for(j=0; j<nvertex; ++j) {
      BZ.push_back(0.0);
    }
    for(j=0; j<Geometry::background_dimension; ++j) {
      x.push_back(-10.0 + 20.0*RND.drandom());
      xnew.push_back(0.0);
    }
    for(j=1+i; j<nvertex; ++j) {
      qt = index.find(make_key(i,j));
      Delta[in1+j-(1+i)] = std::sqrt(distances[qt->second]);
      cdistance.push_back(0.0);
    }
  }

  for(i=0; i<nvertex; ++i) {
    in1 = nvertex*i - i*(i+1)/2;
    for(j=1+i; j<nvertex; ++j) {
      sum = 0.0;
      for(k=0; k<Geometry::background_dimension; ++k) {
        sum += (x[i*Geometry::background_dimension+k]-x[j*Geometry::background_dimension+k])*(x[i*Geometry::background_dimension+k]-x[j*Geometry::background_dimension+k]);
      }
      sum = std::sqrt(sum);
      cdistance[in1+j-(i+1)] = sum;
      lambda = Delta[in1+j-(i+1)];
      err_old += (sum - lambda)*(sum - lambda);
    }
  }

  do {
    for(i=0; i<nvertex; ++i) {
      in1 = nvertex*i - i*(i+1)/2;
      for(j=1+i; j<nvertex; ++j) {
        delta = Delta[in1+j-(i+1)];
        b = -delta/cdistance[in1+j-(i+1)];
        BZ[i*nvertex+j] = b;
        BZ[j*nvertex+i] = b;
      }
    }
    for(i=0; i<nvertex; ++i) {
      sum = 0.0;
      for(j=0; j<nvertex; ++j) {
        if (i == j) continue;
        sum += BZ[nvertex*i+j];
      }
      BZ[nvertex*i+i] = -sum;
    }
    for(i=0; i<nvertex; ++i) {
      for(j=0; j<Geometry::background_dimension; ++j) {
        sum = 0.0;
        for(k=0; k<nvertex; ++k) {
          sum += BZ[i*nvertex+k]*x[Geometry::background_dimension*k+j];
        }
        xnew[Geometry::background_dimension*i+j] = pfactor*sum;
      }
    }
    err = 0.0;
    for(i=0; i<nvertex; ++i) {
      in1 = nvertex*i - i*(i+1)/2;
      for(j=1+i; j<nvertex; ++j) {
        sum = 0.0;
        for(k=0; k<Geometry::background_dimension; ++k) {
          sum += (xnew[i*Geometry::background_dimension+k]-xnew[j*Geometry::background_dimension+k])*(xnew[i*Geometry::background_dimension+k]-xnew[j*Geometry::background_dimension+k]);
        }
        sum = std::sqrt(sum);
        cdistance[in1+j-(i+1)] = sum;
        lambda = Delta[in1+j-(i+1)];
        err += (sum - lambda)*(sum - lambda);
      }
    }
    if (err < cutoff || std::abs(err - err_old) < cutoff) break;
    x = xnew;
    err_old = err;
    its++;
  } while(its <= M);
  */
  int info,nv = nvertex,nwork = 3*nvertex - 1;
  char jtype='V',tsp='N',uplo='U';
  double zero = 0.0,alpha = 1.0;
  std::vector<int> usable;
  hash_map::const_iterator qt;
  const double pfactor = 1.0/double(nvertex);
  double* J = new double[nvertex*nvertex];
  double* A = new double[nvertex*nvertex];
  double* B = new double[nvertex*nvertex];
  double* D = new double[nvertex*nvertex];
  double* work = new double[nwork];
  double* w = new double[nvertex];

  for(i=0; i<nvertex; ++i) {
    usable.push_back(0);
    for(j=0; j<nvertex; ++j) {
      alpha = (j == i) ? 1.0-pfactor : -pfactor;
      J[nvertex*i+j] = alpha;
      qt = index.find(make_key(i,j));
      D[nvertex*i+j] = distances[qt->second];
    }
  }
  // Now form the matrix B = -0.5*J*D2*J
  alpha = 1.0;
  dgemm_(&tsp,&tsp,&nv,&nv,&nv,&alpha,J,&nv,D,&nv,&zero,A,&nv);
  alpha = -0.5;
  dgemm_(&tsp,&tsp,&nv,&nv,&nv,&alpha,A,&nv,J,&nv,&zero,B,&nv);
  // Now compute the eigenvalues and eigenvectors of B
  dsyev_(&jtype,&uplo,&nv,B,&nv,w,work,&nwork,&info);
  if (info != 0) {
    std::cerr << "Error in coordinate calculation, exiting..." << std::endl;
    std::exit(1);
  }
  for(i=0; i<nvertex; ++i) {
    if (w[i] > 0.01) {
      edim++;
      usable[i] = 1;
    }
  }
  std::vector<double> lmatrix,qmatrix;
  for(i=0; i<edim; ++i) {
    for(j=0; j<edim; ++j) {
      lmatrix.push_back(0.0);
    }
  }
  k = 0;
  for(i=0; i<nvertex; ++i) {
    if (w[i] > 0.01) {
      lmatrix[edim*k+k] = std::sqrt(w[i]);
      k++;
    }
  }
  // Transpose the LAPACK output...
  for(i=0; i<nvertex; ++i) {
    for(j=0; j<nvertex; ++j) {
      J[nvertex*i+j] = B[nvertex*j+i];
    }
  }
  for(i=0; i<nvertex; ++i) {
    for(j=0; j<nvertex; ++j) {
      if (usable[j] == 1) qmatrix.push_back(J[nvertex*i+j]);
    }
  }
  for(i=0; i<nvertex; ++j) {
    for(j=0; j<edim; ++j) {
      x.push_back(qmatrix[edim*i+j]*lmatrix[edim*j+j]);
    }
  }
  delete[] A;
  delete[] B;
  delete[] D;
  delete[] J;
  delete[] w;
  delete[] work;

  return edim;
}

double geometry_change(const Geometry* g1,const Geometry* g2)
{
  int i,j,nva;
  bool arg1 = true;
  double gdelta = 0.0;

  assert(g1->relational == g2->relational);
  assert(g1->euclidean == g2->euclidean);

  if (g1->nvertex < g2->nvertex) {
    nva = g1->nvertex;
    arg1 = false;
  }
  else {
    nva = g2->nvertex;
  }
  if (g1->relational) {
    double d1,d2;
    hash_map::const_iterator qt;

    for(i=0; i<nva; ++i) {
      for(j=1+i; j<nva; ++j) {
        qt = g1->index.find(make_key(i,j));
        d1 = g1->distances[qt->second];

        qt = g2->index.find(make_key(i,j));
        d2 = g2->distances[qt->second];
        gdelta += std::abs(d1 - d2);
      }
    }
    if (arg1) {
      for(i=nva; i<g1->nvertex; ++i) {
        for(j=0; j<nva; ++j) {
          qt = g1->index.find(make_key(i,j));
          gdelta += std::abs(g1->distances[qt->second]);
        }
        for(j=1+i; j<g1->nvertex; ++j) {
          qt = g1->index.find(make_key(i,j));
          gdelta += std::abs(g1->distances[qt->second]);
        }
      }
    }
    else {
      for(i=nva; i<g2->nvertex; ++i) {
        for(j=0; j<nva; ++j) {
          qt = g2->index.find(make_key(i,j));
          gdelta += std::abs(g2->distances[qt->second]);
        }
        for(j=1+i; j<g2->nvertex; ++j) {
          qt = g2->index.find(make_key(i,j));
          gdelta += std::abs(g2->distances[qt->second]);
        }
      }
    }
  }
  else {
    double d;
    int n1,n2,m;
    for(i=0; i<nva; ++i) {
      n1 = (signed) g1->coordinates[i].size();
      n2 = (signed) g2->coordinates[i].size();
      m = n2 - n1;
      d = 0.0;
      if (m == 0) {
        for(j=0; j<n1; ++j) {
          d += std::abs(g1->coordinates[i][j] - g2->coordinates[i][j]);
        }
      }
      else if (m > 0) {
        for(j=0; j<n1; ++j) {
          d += std::abs(g1->coordinates[i][j] - g2->coordinates[i][j]);
        }
        for(j=n1; j<n2; ++j) {
          d += std::abs(g2->coordinates[i][j]);
        }
      }
      else {
        for(j=0; j<n2; ++j) {
          d += std::abs(g1->coordinates[i][j] - g2->coordinates[i][j]);
        }
        for(j=n2; j<n1; ++j) {
          d += std::abs(g1->coordinates[i][j]);
        }
      }
      gdelta += d;
    }
    if (arg1) {
      for(i=nva; i<g1->nvertex; ++i) {
        d = 0.0;
        for(j=0; j<(signed) g1->coordinates[i].size(); ++j) {
          d += std::abs(g1->coordinates[i][j]);
        }
        gdelta += d;
      }
    }
    else {
      for(i=nva; i<g2->nvertex; ++i) {
        d = 0.0;
        for(j=0; j<(signed) g2->coordinates[i].size(); ++j) {
          d += std::abs(g2->coordinates[i][j]);
        }
        gdelta += d;
      }
    }
  }
  return gdelta;
}
