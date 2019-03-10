#include "geometry.h"

using namespace SYNARMOSMA;

extern Random RND;

Geometry::Geometry()
{

}

Geometry::Geometry(const Geometry& source)
{
  nvertex = source.nvertex;
  euclidean = source.euclidean;
  relational = source.relational;
  background_dimension = source.background_dimension;
  uniform = source.uniform;
  vperturb = source.vperturb;
  original = source.original;
  distances = source.distances;
  if (relational) {
    coordinates.clear();
  }
  else {
    coordinates = source.coordinates;
  }
}

Geometry& Geometry::operator =(const Geometry& source)
{
  if (this == &source) return *this;

  nvertex = source.nvertex;
  euclidean = source.euclidean;
  relational = source.relational;
  background_dimension = source.background_dimension;
  uniform = source.uniform;
  vperturb = source.vperturb;
  original = source.original;
  distances = source.distances;
  if (relational) {
    coordinates.clear();
  }
  else {
    coordinates = source.coordinates;
  }

  return *this;
}

Geometry::~Geometry()
{

}

void Geometry::initialize(bool type,bool model,bool flat,bool hmemory,int D)
{
  if (D < 1) throw std::invalid_argument("The background dimension of the geometry must be greater than zero!");
  clear();
  euclidean = type;
  relational = model;
  uniform = flat;
  high_memory = hmemory;
  background_dimension = D;
}

void Geometry::clear()
{
  if (!relational) {
    for(int i=0; i<nvertex; ++i) {
      coordinates[i].clear();
    }
    coordinates.clear();
  }
  distances.clear();
  original.clear();
  nvertex = 0;
}

bool Geometry::consistent() const
{
  int i,j;
  const int nc = (signed) coordinates.size();
  const int nd = (signed) distances.size();

  if (nvertex != nc) {
    std::cout << "Illegal length for Geometry::coordinates " << nc << "  " << nvertex << std::endl;
    return false;
  }

  if (!relational) {
    for(i=0; i<nvertex; ++i) {
      for(j=0; j<nc; ++j) {
        if (std::isnan(coordinates[i][j])) {
          std::cout << "NaN in Geometry::coordinates at " << i << " and " << j << std::endl;
          return false;
        }
      }
    }
  }

  for(i=0; i<nvertex; ++i) {
    for(j=1+i; j<nvertex; ++j) {
      if (std::isnan(get_squared_distance(i,j,false))) {
        std::cout << "NaN in Geometry::get_squared_distance for " << i << "  " << j << std::endl;
        std::cout << coordinates[i].size() << "  " << coordinates[j].size() << std::endl;
        return false;
      }
    }
  }

  if (!high_memory && !relational) return true;

  if ((nvertex*(nvertex-1))/2 != nd) {
    std::cout << "Illegal length for Geometry::distances " << nd << "  " << nvertex << "  " << (nvertex*(nvertex-1))/2 << std::endl;
    return false;
  }

  for(i=0; i<nd; ++i) {
    if (std::isnan(distances[i])) {
      std::cout << "NaN in Geometry::distances at " << i << "  " << distances[i] << std::endl;
      return false;
    }
  }
  return true;
}

int Geometry::serialize(std::ofstream& s) const
{
  int i,l,count = 0;
  unsigned int j,n = 0;
#ifdef DISCRETE
  INT64 x;
#else
  double x;
#endif

  s.write((char*)(&nvertex),sizeof(int)); count += sizeof(int);
  s.write((char*)(&background_dimension),sizeof(int)); count += sizeof(int);
  s.write((char*)(&euclidean),sizeof(bool)); count += sizeof(bool);
  s.write((char*)(&relational),sizeof(bool)); count += sizeof(bool);
  s.write((char*)(&uniform),sizeof(bool)); count += sizeof(bool);
  s.write((char*)(&high_memory),sizeof(bool)); count += sizeof(bool);

  if (relational) {
    for(i=0; i<nvertex; ++i) {
      for(l=1+i; l<nvertex; ++l) {
        x = distances[n];
#ifdef DISCRETE
        s.write((char*)(&x),sizeof(INT64)); count += sizeof(INT64);
#else
        s.write((char*)(&x),sizeof(double)); count += sizeof(double);
#endif
        n++;
      }
    }
  }
  else {
    if (uniform) {
      for(i=0; i<nvertex; ++i) {
        for(j=0; j<background_dimension; ++j) {
          x = coordinates[i][j];
#ifdef DISCRETE
          s.write((char*)(&x),sizeof(INT64)); count += sizeof(INT64);
#else
          s.write((char*)(&x),sizeof(double)); count += sizeof(double);
#endif
        }
      }
    }
    else {
      for(i=0; i<nvertex; ++i) {
        n = coordinates[i].size();
        s.write((char*)(&n),sizeof(int)); count += sizeof(int);
        for(j=0; j<n; ++j) {
          x = coordinates[i][j];
#ifdef DISCRETE
          s.write((char*)(&x),sizeof(INT64)); count += sizeof(INT64);
#else
          s.write((char*)(&x),sizeof(double)); count += sizeof(double);
#endif
        }
      }
    }
    if (high_memory) {
      for(i=0; i<nvertex; ++i) {
        for(l=1+i; l<nvertex; ++l) {
          x = distances[n];
#ifdef DISCRETE
          s.write((char*)(&x),sizeof(INT64)); count += sizeof(INT64);
#else
          s.write((char*)(&x),sizeof(double)); count += sizeof(double);
#endif
          n++;
        }
      }
    }
  } 
  return count;
}

int Geometry::deserialize(std::ifstream& s)
{
  int i,j,n,count = 0;
  unsigned int k;
#ifdef DISCRETE
  INT64 x;
#else
  double x;
#endif

  clear();

  s.read((char*)(&nvertex),sizeof(int)); count += sizeof(int);
  s.read((char*)(&background_dimension),sizeof(int)); count += sizeof(int);
  s.read((char*)(&euclidean),sizeof(bool)); count += sizeof(bool);
  s.read((char*)(&relational),sizeof(bool)); count += sizeof(bool);
  s.read((char*)(&uniform),sizeof(bool)); count += sizeof(bool);
  s.read((char*)(&high_memory),sizeof(bool)); count += sizeof(bool);

  if (relational) {
    for(i=0; i<nvertex; ++i) {
      for(j=1+i; j<nvertex; ++j) {
#ifdef DISCRETE
        s.read((char*)(&x),sizeof(INT64)); count += sizeof(INT64);
#else
        s.read((char*)(&x),sizeof(double)); count += sizeof(double);
#endif
        distances.push_back(x);
      }
    }
  }
  else {
#ifdef DISCRETE
    std::vector<INT64> xc;
#else
    std::vector<double> xc;
#endif
    if (uniform) {
      for(i=0; i<nvertex; ++i) {
        for(k=0; k<background_dimension; ++k) {
#ifdef DISCRETE
          s.read((char*)(&x),sizeof(INT64)); count += sizeof(INT64);
#else           
          s.read((char*)(&x),sizeof(double)); count += sizeof(double);
#endif
          xc.push_back(x);
        }
        coordinates.push_back(xc);
        xc.clear();
      }
    }
    else {
      for(i=0; i<nvertex; ++i) {
        s.read((char*)(&n),sizeof(int)); count += sizeof(int);
        for(j=0; j<n; ++j) {
#ifdef DISCRETE
          s.read((char*)(&x),sizeof(INT64)); count += sizeof(INT64);
#else
          s.read((char*)(&x),sizeof(double)); count += sizeof(double);
#endif
          xc.push_back(x);
        }
        coordinates.push_back(xc);
        xc.clear();
      }
    }
    if (high_memory) {
      for(i=0; i<nvertex; ++i) {
        for(j=1+i; j<nvertex; ++j) {
#ifdef DISCRETE
          s.read((char*)(&x),sizeof(INT64)); count += sizeof(INT64);
#else
          s.read((char*)(&x),sizeof(double)); count += sizeof(double);
#endif
          distances.push_back(x);
        }
      }
    }
  }
  return count;
}

double Geometry::perceptual_divergence(const double* raxis,double theta,const double* translation,const double* observed) const
{
  // The method assumes that
  // a) raxis is a unit vector 
  // b) 0 <= theta < 2*pi
  // c) dimension = 3
  if (!uniform || !euclidean || relational || background_dimension != 3) throw std::invalid_argument("The geometry is of the wrong type for the Geometry::perceptual_divergence method!");

  int i,j;
  double ct,st,d1,xt,q0,temp[3],delta,sum = 0.0;

  ct = std::cos(theta);
  st = std::sin(theta);
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i,j,xt,temp,delta,d1,q0) reduction(+:sum)
#endif
  for(i=0; i<nvertex; ++i) {
    // The angular motion:
    // We calculate the cross product...
    temp[0] = coordinates[i][1]*raxis[2] - coordinates[i][2]*raxis[1];
    temp[1] = coordinates[i][2]*raxis[0] - coordinates[i][0]*raxis[0];
    temp[2] = coordinates[i][0]*raxis[1] - coordinates[i][1]*raxis[0];
    // And now the scalar product
    q0 = coordinates[i][0]*raxis[0] + coordinates[i][1]*raxis[1] + coordinates[i][2]*raxis[2];
    delta = 0.0;
    for(j=0; j<3; ++j) {
      xt = q0*raxis[j] + st*temp[j] + ct*(coordinates[i][j] - q0*raxis[j]) + translation[j];
      d1 = xt - observed[3*i+j];
      delta += d1*d1;
    }
    sum += std::sqrt(delta);
  }
  return (sum/double(nvertex));
}

void Geometry::multiple_vertex_addition(int N,bool unf_rnd,const std::vector<double>& x)
{
  if (relational) throw std::invalid_argument("The Geometry::multiple_vertex_addition method is not compatible with a relational geometry!");
  int i,j;
  unsigned int k;
  std::vector<double> xc;
#ifdef DISCRETE
  std::vector<INT64> xi;
#endif
  const unsigned int vsize = x.size();

  clear();

  for(k=0; k<background_dimension; ++k) {
    xc.push_back(0.0);
#ifdef DISCRETE
    xi.push_back(0);
#endif
  }

  if (unf_rnd) {
    if (2*background_dimension != vsize) throw std::invalid_argument("The vector argument in Geometry::multiple_vertex_addition has the wrong length!");

    for(i=0; i<N; ++i) {
      for(k=0; k<background_dimension; ++k) {
        xc[k] = x[2*k] + (x[2*k+1] - x[2*k])*RND.drandom();
      }
#ifdef DISCRETE
      for(k=0; k<background_dimension; ++k) {
        xi[k] = INT64(xc[k]/space_quantum);
      }
      coordinates.push_back(xi);
#else
      coordinates.push_back(xc);
#endif
    }    
  }
  else {
    if (N*background_dimension != vsize) throw std::invalid_argument("The vector argument in Geometry::multiple_vertex_addition has the wrong length!");

    for(i=0; i<N; ++i) {
      for(k=0; k<background_dimension; ++k) {
        xc[k] = x[background_dimension*i + k];
      }
#ifdef DISCRETE
      for(k=0; k<background_dimension; ++k) {
        xi[k] = INT64(xc[k]/space_quantum);
      }
      coordinates.push_back(xi);
#else
      coordinates.push_back(xc);
#endif
    }
  }
  nvertex = N;
  if (!high_memory) return;

  const int npair = N*(N-1)/2;
#ifdef DISCRETE
  INT64 delta;
  const int pfactor = (euclidean) ? 1 : -1;
#else
  double delta;
  const double pfactor = (euclidean) ? 1.0 : -1.0;
#endif

  for(i=0; i<npair; ++i) {
#ifdef DISCRETE
    distances.push_back(0);
#else
    distances.push_back(0.0);
#endif
  }
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i,j,k,delta) schedule(dynamic,1)
#endif
  for(i=0; i<N; ++i) {
    for(j=1+i; j<N; ++j) {
      delta = pfactor*(coordinates[i][0] - coordinates[j][0])*(coordinates[i][0] - coordinates[j][0]);
      for(k=1; k<background_dimension; ++k) {
        delta += (coordinates[i][k] - coordinates[j][k])*(coordinates[i][k] - coordinates[j][k]);
      }
      distances[compute_index(i,j)] = delta;
    }
  }
}

void Geometry::multiple_vertex_addition(int N,double mu,double sigma)
{
  if (relational) throw std::invalid_argument("The Geometry::multiple_vertex_addition method is not compatible with a relational geometry!");
  int i,j;
  unsigned int k;
  std::vector<double> xc;
#ifdef DISCRETE
  std::vector<INT64> xi;
#endif

  clear();

  for(k=0; k<background_dimension; ++k) {
    xc.push_back(0.0);
#ifdef DISCRETE
    xi.push_back(0);
#endif
  }
  for(i=0; i<N; ++i) {
    for(k=0; k<background_dimension; ++k) {
      xc[k] = RND.nrandom(mu,sigma);
    }
#ifdef DISCRETE
    for(k=0; k<background_dimension; ++k) {
      xi[k] = INT64(xc[k]/space_quantum);
    }
    coordinates.push_back(xi);
#else
    coordinates.push_back(xc);
#endif
  }
  nvertex = N;
  if (!high_memory) return;

  const int npair = N*(N-1)/2;
#ifdef DISCRETE
  INT64 delta;
  const int pfactor = (euclidean) ? 1 : -1;
#else
  double delta;
  const double pfactor = (euclidean) ? 1.0 : -1.0;
#endif

  for(i=0; i<npair; ++i) {
#ifdef DISCRETE
    distances.push_back(0);
#else
    distances.push_back(0.0);
#endif
  }

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i,j,k,delta) schedule(dynamic,1)
#endif
  for(i=0; i<N; ++i) {
    for(j=1+i; j<N; ++j) {
      delta = pfactor*(coordinates[i][0] - coordinates[j][0])*(coordinates[i][0] - coordinates[j][0]);
      for(k=1; k<background_dimension; ++k) {
        delta += (coordinates[i][k] - coordinates[j][k])*(coordinates[i][k] - coordinates[j][k]);
      }
      distances[compute_index(i,j)] = delta;
    }
  }
}

void Geometry::multiple_vertex_addition(const std::vector<std::vector<double> >& source)
{
  if (relational) throw std::invalid_argument("The Geometry::multiple_vertex_addition method is not compatible with a relational geometry!");
  int i,j;
  unsigned int k;

  clear();

  nvertex = (signed) source.size();
#ifdef DISCRETE
  INT64 n;
  std::vector<INT64> xi;
  for(i=0; i<nvertex; ++i) {
    for(k=0; k<source[i].size(); ++k) {
      n = INT64(source[i][k]/space_quantum);
      xi.push_back(n);
    }
    coordinates.push_back(xi);
    xi.clear();
  }
#else
  coordinates = source;
#endif

  if (!high_memory) return;

  const int npair = nvertex*(nvertex-1)/2;
#ifdef DISCRETE
  INT64 delta;
  const int pfactor = (euclidean) ? 1 : -1;
#else
  double delta;
  const double pfactor = (euclidean) ? 1.0 : -1.0;
#endif

  for(i=0; i<npair; ++i) {
#ifdef DISCRETE
    distances.push_back(0);
#else
    distances.push_back(0.0);
#endif
  }
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i,j,k,delta) schedule(dynamic,1)
#endif
  for(i=0; i<nvertex; ++i) {
    for(j=1+i; j<nvertex; ++j) {
      delta = pfactor*(coordinates[i][0] - coordinates[j][0])*(coordinates[i][0] - coordinates[j][0]);
      for(k=1; k<background_dimension; ++k) {
        delta += (coordinates[i][k] - coordinates[j][k])*(coordinates[i][k] - coordinates[j][k]);
      }
      distances[compute_index(i,j)] = delta;
    }
  }
}

double Geometry::dot_product(const std::vector<double>& vx,const std::vector<double>& vy) const
{
  double output = 0.0;
  unsigned int i;
  if (uniform) {
    output = (euclidean) ? vx[0]*vy[0] : -vx[0]*vy[0];
    for(i=1; i<background_dimension; ++i) {
      output += vx[i]*vy[i];
    }
    return output;
  }
  unsigned int l,k = 0,n = vx.size(),m = vy.size();
  std::vector<double> vlx = vx;
  std::vector<double> vly = vy;

  l = n;

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
  output = (euclidean) ? vx[0]*vy[0] : -vx[0]*vy[0];
  for(i=0; i<l; ++i) {
    output += vlx[i]*vly[i];
  }
  return output;
}

void Geometry::get_implied_vertices(int n,std::set<int>& vx) const
{
  int i,j;
  vx.clear();
  if (relational) {
    // Loop over all vertex pairs and find the one that corresponds to index
    // "n"
    for(i=0; i<nvertex; ++i) {
      for(j=1+i; j<nvertex; ++j) {
        if (compute_index(i,j) == n) {
          vx.insert(i);
          vx.insert(j);
          return;
        }
      }
    }
  }
  else {
    if (uniform) {
      vx.insert(n/background_dimension);
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
}

void Geometry::load(const Geometry* source)
{
  nvertex = source->nvertex;
  euclidean = source->euclidean;
  relational = source->relational;
  uniform = source->uniform;
  vperturb = source->vperturb;
  background_dimension = source->background_dimension;
  distances = source->distances;
  if (!relational) coordinates = source->coordinates;
}

void Geometry::store(Geometry* target) const
{
  target->nvertex = nvertex;
  target->euclidean = euclidean;
  target->relational = relational;
  target->uniform = uniform;
  target->vperturb = vperturb;
  target->background_dimension = background_dimension;
  target->distances = distances;
  if (!relational) target->coordinates = coordinates;
}

void Geometry::create(int n,std::string& type)
{
  if (!relational) throw std::runtime_error("Illegal geometric method (initialize) call for absolute model!");

  distances.clear();

  int i,j;
  boost::to_upper(type);
  if (type == "CARTESIAN") {
    nvertex = ipow(n,background_dimension);
    if (!high_memory) return;

    int p,q,s,alpha;
    unsigned int l;
#ifdef DISCRETE
    INT64 delta;
#else
    double delta;
#endif
    std::vector<int> x,y;

    for(i=0; i<nvertex; ++i) {
      s = i;
      for(l=background_dimension; l>=1; l--) {
        p = ipow(n,l-1);
        q = s/p;
        x.push_back(q);
        s -= p*q;
      }
      x.push_back(s);
      for(j=1+i; j<nvertex; ++j) {
        s = j;
        for(l=background_dimension; l>=1; --l) {
          p = ipow(n,l-1);
          q = s/p;
          y.push_back(q);
          s -= p*q;
        }
        y.push_back(s);
#ifdef DISCRETE
        delta = INT64((x[0] - y[0])*(x[0] - y[0]));
#else
        delta = double((x[0] - y[0])*(x[0] - y[0]));
#endif
        if (!euclidean) delta = -delta;
        for(l=1; l<background_dimension; ++l) {
          alpha = (x[l] - y[l])*(x[l] - y[l]);
#ifdef DISCRETE
          delta += INT64(alpha);
#else
          delta += double(alpha);
#endif
        }
        distances.push_back(delta);
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
    if (!high_memory) return;
    if (euclidean) {
      for(i=0; i<nvertex; ++i) {
        for(j=1+i; j<nvertex; ++j) {
#ifdef DISCRETE
          distances.push_back(1 + INT64(RND.irandom(250)));
#else
          distances.push_back(0.25 + 10.0*RND.drandom());
#endif
        }
      }
    }
    else {
#ifdef DISCRETE
      INT64 alpha;
#else
      double alpha,r = 1.0 + double(background_dimension - 1);
#endif
      for(i=0; i<nvertex; ++i) {
        for(j=1+i; j<nvertex; ++j) {
#ifdef DISCRETE
          alpha = -5 + 5*INT64(RND.irandom(background_dimension - 1));
#else
          alpha = -5.0 + 5.0*r*RND.drandom();
#endif
          distances.push_back(alpha);
        }
      }
    }
  }
  else if (type == "MONOPLEX") {
    // In this case n is the monoplex dimension, so the number of vertices
    // is 1+n
    nvertex = 1 + n;
    if (!high_memory) return;
#ifdef DISCRETE
    const INT64 zero = 0;
    const INT64 one = 1;
    const INT64 m_one = -1;
    const INT64 two = 2;
#else
    const double zero = 0.0;
    const double one = 1.0;
    const double m_one = -1.0;
    const double two = 2.0;
#endif
    if (euclidean) {
      // The distance from vertex 0 to the other vertices is unity...
      for(i=0; i<nvertex-1; ++i) {
        distances.push_back(one);
      }
      // All the other distances are equal to 2.0
      for(i=1; i<nvertex; ++i) {
        for(j=1+i; j<nvertex; ++j) {
          distances.push_back(two);
        }
      }
    }
    else {
      // Distances between v_0 and the other n-1 vertices
      distances.push_back(m_one);
      for(i=0; i<nvertex-2; ++i) {
        distances.push_back(one);
      }
      // Distances between v_1 and the other n-2 vertices
      for(i=0; i<nvertex-2; ++i) {
        distances.push_back(zero);
      }
      // All the other distances...
      for(i=2; i<nvertex; ++i) {
        for(j=1+i; j<nvertex; ++j) {
          distances.push_back(two);
        }
      }
    }
  }
  else {
    throw std::invalid_argument("Unrecognized spatial type in Geometry::create method!");
  }
}

void Geometry::vertex_difference(int n,int m,std::vector<double>& delta) const
{
  if (n == m) throw std::invalid_argument("The vertex arguments in Geometry::vertex_difference must be distinct!");
  if (relational) throw std::runtime_error("Illegal geometric method (vertex_difference) call for relational model!");

  unsigned int i;

  delta.clear();
  if (uniform) {
    double xd;
    for(i=0; i<background_dimension; ++i) {
#ifdef DISCRETE
      xd = space_quantum*double(coordinates[n][i] - coordinates[m][i]);
#else
      xd = coordinates[n][i] - coordinates[m][i];
#endif
      delta.push_back(xd);
    }
    return;
  }
  unsigned int d1 = coordinates[n].size();
  unsigned int d2 = coordinates[m].size();
  unsigned int k = 0,l = d1;
  std::vector<double> vlx,vly;

#ifdef DISCRETE
  for(i=0; i<d1; ++i) {
    vlx.push_back(space_quantum*double(coordinates[n][i]));
  }
  for(i=0; i<d2; ++i) {
    vly.push_back(space_quantum*double(coordinates[m][i]));
  }
#else
  vlx = coordinates[n];
  vly = coordinates[m];
#endif

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

void Geometry::rollback(bool minimal)
{
  if (vperturb == -1) std::runtime_error("An unperturbed geometry cannot be rolled back!");

  if (minimal) {
    if (relational) {
      distances[vperturb] = original[0];
    }
    else {
      int i,kt = 0;
      unsigned int j;

      for(i=0; i<nvertex; ++i) {
        for(j=0; j<coordinates[i].size(); ++j) {
          if (kt == vperturb) {
            coordinates[i][j] = original[0];
            return;
          }
          kt++;
        }
      }
    }
  }
  else {
    if (relational) {
      int i,j = 0;
      for(i=0; i<nvertex; ++i) {
        if (i == vperturb) continue;
        distances[compute_index(i,vperturb)] = original[j];
        j++;
      }
    }
    else {
      coordinates[vperturb] = original;
    }
  }
  vperturb = -1;
  original.clear();
}

int Geometry::vertex_addition(const std::set<int>& antecedents)
{
  if (antecedents.empty()) {
    vertex_addition(-1);
  }
  else {
    if (RND.drandom() < 0.25) {
      return vertex_addition(RND.irandom(antecedents));
    }
    else {
      int i,j;
      std::set<int>::const_iterator it;
      const double na = double(antecedents.size());

      if (relational) {
        int k = 0;
#ifdef DISCRETE
        INT64 l;
        std::vector<INT64> ndistances;
        const INT64 nt = INT64(antecedents.size());
#else
        double l;
        std::vector<double> ndistances;
#endif

        for(i=0; i<nvertex; ++i) {
          for(j=1+i; j<nvertex; ++j) {
            ndistances.push_back(distances[k]);
            k++;
          }
#ifdef DISCRETE
          l = 0;
          for(it=antecedents.begin(); it!=antecedents.end(); ++it) {
            l += distances[compute_index(i,*it)];
          }
          l = l/nt;
#else
          l = 0.0;
          for(it=antecedents.begin(); it!=antecedents.end(); ++it) {
            l += distances[compute_index(i,*it)];
          }
          l = l/na;
#endif
          ndistances.push_back(l);
        }
        distances = ndistances;
        nvertex++;
        return nvertex;
      }
      else {
        unsigned int k;
        std::vector<double> xc,avg_x;
        for(k=0; k<background_dimension; ++k) {
          avg_x.push_back(0.0);
        }
        for(it=antecedents.begin(); it!=antecedents.end(); ++it) {
          j = *it;
          for(k=0; k<background_dimension; ++k) {
#ifdef DISCRETE
            avg_x[k] += space_quantum*double(coordinates[j][k]);
#else
            avg_x[k] += coordinates[j][k];
#endif
          }
        }
        for(k=0; k<background_dimension; ++k) {
          avg_x[k] /= na;
        }
        for(k=0; k<background_dimension; ++k) {
          xc.push_back(RND.nrandom(avg_x[k],0.5));
        }
        return vertex_addition(xc);
      }
    }
  }
  return nvertex;
}

int Geometry::vertex_addition(int parent,double mutation)
{
  int i,j,k = 0;
  double alpha;
  std::vector<double> x;

  if (parent == -1) {
    // No antecedent, so place the vertex randomly...
    if (relational) {
#ifdef DISCRETE
      std::vector<INT64> ndistances;
#else
      std::vector<double> ndistances;
#endif

      if (euclidean) {
        for(i=0; i<nvertex; ++i) {
          for(j=1+i; j<nvertex; ++j) {
            ndistances.push_back(distances[k]);
            k++;
          }
          alpha = 1.0 + 15.0*RND.drandom();
          ndistances.push_back(alpha);
        }
      }
      else {
        // Skew the random number generation so the inter-vertex
        // lengths are properly weighted between timelike and
        // spacelike...
        double r = 1.0 + double(background_dimension-1);
        for(i=0; i<nvertex; ++i) {
          for(j=1+i; j<nvertex; ++j) {
            ndistances.push_back(distances[k]);
            k++;
          }
          alpha = -5.0 + 5.0*r*RND.drandom();
          ndistances.push_back(alpha);
        }
      }
      distances = ndistances;
      nvertex++;
    }
    else {
      unsigned int l;

      for(l=0; l<background_dimension; ++l) {
        alpha = -10.0 + 20.0*RND.drandom();
        x.push_back(alpha);
      }
      vertex_addition(x);
    }
  }
  else {
    // Should be close to its parent vertex...
    if (relational) {
      double mu,sigma;
#ifdef DISCRETE
      std::vector<INT64> ndistances;
#else
      std::vector<double> ndistances;
#endif

      for(i=0; i<nvertex; ++i) {
        for(j=1+i; j<nvertex; ++j) {
          ndistances.push_back(distances[k]);
          k++;
        }
        if (i == parent) {
          mu = 0.0; 
          sigma = 0.1;
        }
        else {
          mu = distances[compute_index(i,parent)];
          sigma = mutation/10.0;
        }
        ndistances.push_back(RND.nrandom(mu,sigma));
      }
      distances = ndistances;
      nvertex++;
    }
    else {
      unsigned int q,p = RND.irandom(background_dimension);
      double r = RND.drandom(0.1+0.5*mutation,0.2+mutation);
      get_coordinates(parent,x);
      alpha = RND.drandom(0.0,2.0*M_PI);
      do {
        q = RND.irandom(background_dimension);
        if (q != p) break;
      } while(true);
      x[p] += r*std::cos(alpha);
      x[q] += r*std::sin(alpha);
      vertex_addition(x);
    }
  }
  return nvertex;
}

void Geometry::mutation(int v,bool by_vertex,bool complete,double severity)
{
#ifdef DISCRETE
  INT64 alpha = INT64((RND.nrandom(0.0,severity)/space_quantum));
#else
  double alpha = RND.nrandom(0.0,severity);
#endif

  if (by_vertex) {
    if (v < 0 || v >= nvertex) throw std::invalid_argument("Illegal vertex value in Geometry::mutation method!");
    vperturb = v;
    if (relational) {
      int i,j;

      original.clear();
      for(i=0; i<nvertex; ++i) {
        if (i == v) continue;
        j = compute_index(v,i);
        original.push_back(distances[j]);
#ifdef DISCRETE
        alpha = INT64((RND.nrandom(0.0,severity)/space_quantum));
#else
        alpha = RND.nrandom(0.0,severity);
#endif
        distances[j] += alpha;
      }
    }
    else {
      unsigned int i;

      original = coordinates[v];
      if (complete) {
        for(i=0; i<coordinates[v].size(); ++i) {
#ifdef DISCRETE
          alpha = INT64((RND.nrandom(0.0,severity)/space_quantum));
#else
          alpha = RND.nrandom(0.0,severity);
#endif
          coordinates[v][i] += alpha;
        }
      }
      else {
        i = coordinates[v].size();
        coordinates[v][RND.irandom(i)] += RND.nrandom(0.0,severity);
      }
    }
  }
  else {
    vperturb = v;
    original.clear();
    if (relational) {
      original.push_back(distances[v]);
      distances[v] += alpha;
      return;
    }
    int i,kt = 0;
    unsigned int j;
    for(i=0; i<nvertex; ++i) {
      for(j=0; j<coordinates[i].size(); ++j) {
        if (kt == v) {
          original.push_back(coordinates[i][j]);
          coordinates[i][j] += alpha;
          return;
        }
        kt++;
      }
    }
  }
}

bool Geometry::adjust_dimension(const std::vector<int>& vdimension)
{
  if (nvertex != (signed) vdimension.size()) throw std::invalid_argument("The length of the dimension vector must be equal to the number of vertices!");
  if (relational || uniform) return false;
  
  int i;
  unsigned int j,n,m;
  std::set<int> vmodified;
#ifdef DISCRETE
  std::vector<INT64> x;
#else
  std::vector<double> x;
#endif

  for(i=0; i<nvertex; ++i) {
    if (vdimension[i] == -1) continue;
    // First check to see if the dimension has changed...
    m = (signed) vdimension[i];
    n = coordinates[i].size();
    if (m <= background_dimension) {
      if (n != background_dimension) {
        vmodified.insert(i);
        x = coordinates[i];
        coordinates[i].clear();
        for(j=0; j<background_dimension; ++j) {
          coordinates[i].push_back(x[j]);
        }
      }
    }
    else {
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
  if (!vmodified.empty()) {
    compute_squared_distances(vmodified);
    return true;
  }
  return false;
}

void Geometry::compute_relational_matrices(std::vector<double>& R,std::vector<std::vector<double> >& angles) const
{
  R.clear();
  angles.clear();
  // This method supposes that we are working in Euclidean n-space (n > 1) with a 
  // uniform, absolute geometry
  if (background_dimension < 2 || !euclidean || !uniform || relational) throw std::invalid_argument("The geometry is of the wrong type for the Geometry::compute_relational_matrices method!");

  int i,j,k;
  unsigned int l;
  double r,v[background_dimension];
#ifdef DISCRETE
  INT64 base[background_dimension];
#else
  double base[background_dimension];
#endif

  // For the angles matrices, the last one is always the one that covers the range [0,2\pi), so 
  // all the others range over [0,\pi]
  if (background_dimension == 2) {
    for(i=0; i<nvertex*nvertex; ++i) {
      R.push_back(0.0);
    }
    angles.push_back(R);
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i,j,base,v,r)
#endif
    for(i=0; i<nvertex; ++i) {
      base[0] = coordinates[i][0];
      base[1] = coordinates[i][1];
      for(j=0; j<nvertex; ++j) {
        if (i == j) continue;
#ifdef DISCRETE
        v[0] = space_quantum*double(coordinates[j][0] - base[0]);
        v[1] = space_quantum*double(coordinates[j][1] - base[1]);
#else
        v[0] = coordinates[j][0] - base[0];
        v[1] = coordinates[j][1] - base[1];
#endif
        R[j+nvertex*i] = std::sqrt(get_squared_distance(i,j,false));
        angles[0][j+nvertex*i] = M_PI + std::atan2(v[1],v[0]);
      }
    }
  }
  else if (background_dimension == 3) {
    for(i=0; i<nvertex*nvertex; ++i) {
      R.push_back(0.0);
    }
    angles.push_back(R);
    angles.push_back(R);
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i,j,base,v,r)
#endif
    for(i=0; i<nvertex; ++i) {
      base[0] = coordinates[i][0];
      base[1] = coordinates[i][1];
      base[2] = coordinates[i][2];
      for(j=0; j<nvertex; ++j) {
        if (i == j) continue;
#ifdef DISCRETE
        v[0] = space_quantum*double(coordinates[j][0] - base[0]);
        v[1] = space_quantum*double(coordinates[j][1] - base[1]);
        v[2] = space_quantum*double(coordinates[j][2] - base[2]);
#else
        v[0] = coordinates[j][0] - base[0];
        v[1] = coordinates[j][1] - base[1];
        v[2] = coordinates[j][2] - base[2];
#endif
        r = std::sqrt(get_squared_distance(i,j,false));
        R[j+nvertex*i] = r;
        angles[0][j+nvertex*i] = std::acos(v[2]/r);
        angles[1][j+nvertex*i] = M_PI + std::atan2(v[1],v[0]); 
      }
    }
  }
  else {
    double sum;
    const int nm1 = background_dimension - 1;
    const int nm2 = background_dimension - 2;
    for(i=0; i<nvertex*nvertex; ++i) {
      R.push_back(0.0);
    }
    for(i=0; i<nm1; ++i) {
      angles.push_back(R);
    }
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i,j,k,base,v,r,sum)
#endif
    for(i=0; i<nvertex; ++i) {
      for(l=0; l<background_dimension; ++l) {
        base[l] = coordinates[i][l];
      }
      for(j=0; j<nvertex; ++j) {
        if (i == j) continue;
        for(l=0; l<background_dimension; ++l) {
#ifdef DISCRETE
          v[l] = space_quantum*double(coordinates[j][l] - base[l]);
#else
          v[l] = coordinates[j][l] - base[l];
#endif
        }
        r = get_squared_distance(i,j,false);
        R[j+nvertex*i] = std::sqrt(r);
        sum = 0.0;
        for(k=0; k<nm2; ++k) {
          angles[k][j+nvertex*i] = std::acos(v[k]/std::sqrt(r - sum));
          sum += v[k]*v[k];
        }
        r = std::sqrt(v[nm1]*v[nm1] + v[nm2]*v[nm2]);
        if (v[nm1] < 0.0) {
          angles[nm2][j+nvertex*i] = 2.0*M_PI - std::acos(v[nm2]/r);
        }
        else {
          angles[nm2][j+nvertex*i] = std::acos(v[nm2]/r);
        }
      }
    }
  }
}

void Geometry::compute_squared_distances(const std::set<int>& vmodified)
{
  if (relational || !high_memory) return;

  int i,j,in1;
  unsigned int n1,n2,k;
  std::set<int>::const_iterator it;
#ifdef DISCRETE
  INT64 delta;
  const int pfactor = (euclidean) ? 1 : -1;
#else
  double delta;
  const double pfactor = (euclidean) ? 1.0 : -1.0;
#endif

  if (uniform) {
    for(it=vmodified.begin(); it!=vmodified.end(); ++it) {
      i = *it;
      for(j=0; j<i; ++j) {
        in1 = j*nvertex - j*(j+1)/2;
        delta = pfactor*(coordinates[i][0] - coordinates[j][0])*(coordinates[i][0] - coordinates[j][0]);
        for(k=1; k<background_dimension; ++k) {
          delta += (coordinates[i][k] - coordinates[j][k])*(coordinates[i][k] - coordinates[j][k]);
        }
        distances[in1+i-(1+j)] = delta;
      }
      in1 = i*nvertex - i*(i+1)/2;
      for(j=1+i; j<nvertex; ++j) {
        delta = pfactor*(coordinates[i][0] - coordinates[j][0])*(coordinates[i][0] - coordinates[j][0]);
        for(k=1; k<background_dimension; ++k) {
          delta += (coordinates[i][k] - coordinates[j][k])*(coordinates[i][k] - coordinates[j][k]);
        }
        distances[in1+j-(1+i)] = delta;
      }
    }
  }
  else {
    for(it=vmodified.begin(); it!=vmodified.end(); ++it) {
      i = *it;
      n1 = coordinates[i].size();
      for(j=0; j<i; ++j) {
        in1 = j*nvertex - j*(j+1)/2;
        delta = pfactor*(coordinates[i][0] - coordinates[j][0])*(coordinates[i][0] - coordinates[j][0]);
        n2 = (signed) coordinates[j].size();
        n2 = (n1 <= n2) ? n1 : n2;
        for(k=1; k<n2; ++k) {
          delta += (coordinates[i][k] - coordinates[j][k])*(coordinates[i][k] - coordinates[j][k]);
        }
        distances[in1+i-(1+j)] = delta;
      }
      in1 = i*nvertex - i*(i+1)/2;
      for(j=1+i; j<nvertex; ++j) {
        n2 = coordinates[j].size();
        n2 = (n1 <= n2) ? n1 : n2;
        delta = pfactor*(coordinates[i][0] - coordinates[j][0])*(coordinates[i][0] - coordinates[j][0]);
        for(k=1; k<n2; ++k) {
          delta += (coordinates[i][k] - coordinates[j][k])*(coordinates[i][k] - coordinates[j][k]);
        }
        distances[in1+j-(1+i)] = delta;
      }
    }
  }
}

void Geometry::compute_squared_distances()
{
  if (relational || !high_memory) return;

  int i,j,in1;
  unsigned int k;
  const int npair = nvertex*(nvertex-1)/2;
#ifdef DISCRETE
  INT64 delta;
  const int pfactor = (euclidean) ? 1 : -1;
#else
  double delta;
  const double pfactor = (euclidean) ? 1.0 : -1.0;
#endif

  distances.clear();
  for(i=0; i<npair; ++i) {
#ifdef DISCRETE
    distances.push_back(0);
#else
    distances.push_back(0.0);
#endif
  }

  if (uniform) {
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i,j,k,in1,delta) schedule(dynamic,1)
#endif
    for(i=0; i<nvertex; ++i) {
      in1 = i*nvertex - i*(i+1)/2;
      for(j=1+i; j<nvertex; ++j) {
        delta = pfactor*(coordinates[i][0] - coordinates[j][0])*(coordinates[i][0] - coordinates[j][0]);
        for(k=1; k<background_dimension; ++k) {
          delta += (coordinates[i][k] - coordinates[j][k])*(coordinates[i][k] - coordinates[j][k]);
        }
        distances[in1+j-(1+i)] = delta;
      }
    }
  }
  else {
    unsigned int n1,n2;
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i,j,k,n1,n2,in1,delta) schedule(dynamic,1)
#endif
    for(i=0; i<nvertex; ++i) {
      n1 = coordinates[i].size();
      in1 = i*nvertex - i*(i+1)/2;
      for(j=1+i; j<nvertex; ++j) {
        delta = pfactor*(coordinates[i][0] - coordinates[j][0])*(coordinates[i][0] - coordinates[j][0]);
        n2 = coordinates[j].size();
        n2 = (n1 <= n2) ? n1 : n2;
        for(k=1; k<n2; ++k) {
          delta += (coordinates[i][k] - coordinates[j][k])*(coordinates[i][k] - coordinates[j][k]);
        }
        distances[in1+j-(1+i)] = delta;
      }
    }
  }
}

unsigned int Geometry::compute_coordinates(std::vector<double>& x) const
{
  int i;
  unsigned int k,edim = 0;

  x.clear();

  if (!relational) {
    unsigned int j;
    // Calculate the maximum simplicial dimension...
    for(i=0; i<nvertex; ++i) {
      j = coordinates[i].size();
      if (j > edim) edim = j;
    }
    // Pad with zeroes as necessary...
    for(i=0; i<nvertex; ++i) {
      k = coordinates[i].size();
      for(j=0; j<k; ++j) {
#ifdef DISCRETE 
        x.push_back(space_quantum*double(coordinates[i][j]));
#else
        x.push_back(coordinates[i][j]);
#endif
      }
      for(j=k; j<edim; ++j) {
        x.push_back(0.0);
      }
    }
    return edim;
  }
  if (!euclidean) throw std::runtime_error("If the geometry is relational it must also be Euclidean in the Geometry::compute_coordinates method!");

  if (nvertex < 2) throw std::invalid_argument("The number of vertices must be greater than unity in the Geometry::compute_coordinates method!");

  int j,info,nv = nvertex,nwork = 3*nvertex - 1;
  char jtype='V',tsp='N',uplo='U';
  double zero = 0.0,alpha = 1.0;
  bool usable[nvertex];
  const double pfactor = 1.0/double(nvertex);
  double* J = new double[nvertex*nvertex];
  double* A = new double[nvertex*nvertex];
  double* B = new double[nvertex*nvertex];
  double* D = new double[nvertex*nvertex];
  double* work = new double[nwork];
  double* w = new double[nvertex];

  // Put matrices in column-major form...
  for(i=0; i<nvertex; ++i) {
    usable[i] = false;
    for(j=0; j<nvertex; ++j) {
      alpha = (j == i) ? 1.0 - pfactor : -pfactor;
      J[nvertex*j+i] = alpha;
      D[nvertex*j+i] = std::sqrt(get_squared_distance(i,j,false)); 
    }
  }
  // Now form the matrix B = -0.5*J*D2*J
  alpha = 1.0;
  dgemm_(&tsp,&tsp,&nv,&nv,&nv,&alpha,J,&nv,D,&nv,&zero,A,&nv);
  alpha = -0.5;
  dgemm_(&tsp,&tsp,&nv,&nv,&nv,&alpha,A,&nv,J,&nv,&zero,B,&nv);
  // Now compute the eigenvalues and eigenvectors of B
  dsyev_(&jtype,&uplo,&nv,B,&nv,w,work,&nwork,&info);
  if (info != 0) throw std::runtime_error("Error in Geometry::compute_coordinates calculation!");
  for(i=0; i<nvertex; ++i) {
    if (w[i] > std::numeric_limits<double>::epsilon()) {
      edim++;
      usable[i] = true;
    }
  }
  std::vector<double> lmatrix,qmatrix;
  for(k=0; k<edim*edim; ++k) {
    lmatrix.push_back(0.0);
  }
  k = 0;
  for(i=0; i<nvertex; ++i) {
    if (w[i] > std::numeric_limits<double>::epsilon()) {
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
      if (usable[j]) qmatrix.push_back(J[nvertex*i+j]);
    }
  }
  for(i=0; i<nvertex; ++j) {
    for(k=0; k<edim; ++k) {
      x.push_back(qmatrix[edim*i+k]*lmatrix[edim*k+k]);
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

double SYNARMOSMA::geometry_change(const Geometry* g1,const Geometry* g2)
{
  int i,j,nva;
  bool arg1 = true;
#ifdef DISCRETE
  INT64 gdelta = 0;
#else
  double gdelta = 0.0;
#endif
  if (g1->euclidean != g2->euclidean) throw std::invalid_argument("The two geometries in geometry_change are inconsistent!");
  if (g1->relational != g2->relational) throw std::invalid_argument("The two geometries in geometry_change are inconsistent!");

  if (g1->nvertex < g2->nvertex) {
    nva = g1->nvertex;
    arg1 = false;
  }
  else {
    nva = g2->nvertex;
  }
  if (g1->relational) {
#ifdef DISCRETE
    INT64 d1,d2;
#else
    double d1,d2;
#endif

    for(i=0; i<nva; ++i) {
      for(j=1+i; j<nva; ++j) {
        d1 = g1->distances[g1->compute_index(i,j)];

        d2 = g2->distances[g2->compute_index(i,j)];
        gdelta += std::abs(d1 - d2);
      }
    }
    if (arg1) {
      for(i=nva; i<g1->nvertex; ++i) {
        for(j=0; j<nva; ++j) {
          gdelta += std::abs(g1->distances[g1->compute_index(i,j)]);
        }
        for(j=1+i; j<g1->nvertex; ++j) {
          gdelta += std::abs(g1->distances[g1->compute_index(i,j)]);
        }
      }
    }
    else {
      for(i=nva; i<g2->nvertex; ++i) {
        for(j=0; j<nva; ++j) {
          gdelta += std::abs(g2->distances[g2->compute_index(i,j)]);
        }
        for(j=1+i; j<g2->nvertex; ++j) {
          gdelta += std::abs(g2->distances[g2->compute_index(i,j)]);
        }
      }
    }
  }
  else {
    int n1,n2,m;
#ifdef DISCRETE
    INT64 d;
    const INT64 zero = 0;
#else
    double d;
    const double zero = 0.0;
#endif

    for(i=0; i<nva; ++i) {
      n1 = (signed) g1->coordinates[i].size();
      n2 = (signed) g2->coordinates[i].size();
      m = n2 - n1;
      d = zero;
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
        d = zero;
        for(j=0; j<(signed) g1->coordinates[i].size(); ++j) {
          d += std::abs(g1->coordinates[i][j]);
        }
        gdelta += d;
      }
    }
    else {
      for(i=nva; i<g2->nvertex; ++i) {
        d = zero;
        for(j=0; j<(signed) g2->coordinates[i].size(); ++j) {
          d += std::abs(g2->coordinates[i][j]);
        }
        gdelta += d;
      }
    }
  }
#ifdef DISCRETE
  return space_quantum*double(gdelta);
#else
  return gdelta;
#endif
}
