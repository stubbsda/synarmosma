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

#include "geometry.h"

using namespace SYNARMOSMA;

extern Random RND;

Geometry::Geometry()
{
  set_default_values();
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
  clear();
}

void Geometry::set_default_values()
{
  vperturb = -1;
  nvertex = 0;
  euclidean = true;
  relational = false;
  uniform = true;
  high_memory = true;
  background_dimension = 3;
}

void Geometry::initialize(bool type,bool model,bool flat,bool hmemory,int D)
{
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
#ifdef DEBUG
  assert(nvertex == (signed) coordinates.size());
#endif
  int i,j;

  if (!relational) {
    for(i=0; i<nvertex; ++i) {
      for(j=0; j<(signed) coordinates[i].size(); ++j) {
        if (std::isnan(coordinates[i][j])) {
          std::cout << "Nan at " << i << " and " << j << std::endl;
          return false;
        }
      }
    }
  }

  for(i=0; i<nvertex; ++i) {
    for(j=1+i; j<nvertex; ++j) {
      if (std::isnan(get_distance(i,j,false))) {
        std::cout << i << "  " << j << std::endl;
        std::cout << coordinates[i].size() << "  " << coordinates[j].size() << std::endl;
        return false;
      }
    }
  }

  if (!high_memory) return true;

  if ((nvertex*(nvertex-1))/2 != (signed) distances.size()) {
    std::cout << distances.size() << "  " << nvertex << "  " << (nvertex*(nvertex-1))/2 << std::endl;
    return false;
  }

  for(i=0; i<(signed) distances.size(); ++i) {
    if (std::isnan(distances[i])) {
      std::cout << i << "  " << distances[i] << std::endl;
      return false;
    }
  }
  return true;
}

void Geometry::serialize(std::ofstream& s) const
{
  int i,j,n = 0;
  double x;

  s.write((char*)(&nvertex),sizeof(int));
  s.write((char*)(&background_dimension),sizeof(int));
  s.write((char*)(&euclidean),sizeof(bool));
  s.write((char*)(&relational),sizeof(bool));
  s.write((char*)(&uniform),sizeof(bool));
  s.write((char*)(&high_memory),sizeof(bool));

  if (relational) {
    for(i=0; i<nvertex; ++i) {
      for(j=1+i; j<nvertex; ++j) {
        x = distances[n];
        s.write((char*)(&x),sizeof(double));
        n++;
      }
    }
  }
  else {
    if (uniform) {
      for(i=0; i<nvertex; ++i) {
        for(j=0; j<background_dimension; ++j) {
          x = coordinates[i][j];
          s.write((char*)(&x),sizeof(double));
        }
      }
    }
    else {
      for(i=0; i<nvertex; ++i) {
        n = (signed) coordinates[i].size();
        s.write((char*)(&n),sizeof(int));
        for(j=0; j<n; ++j) {
          x = coordinates[i][j];
          s.write((char*)(&x),sizeof(double));
        }
      }
    }
    if (high_memory) {
      for(i=0; i<nvertex; ++i) {
        for(j=1+i; j<nvertex; ++j) {
          x = distances[n];
          s.write((char*)(&x),sizeof(double));
          n++;
        }
      }
    }
  } 
}

void Geometry::deserialize(std::ifstream& s)
{
  int i,j,n;
#ifdef DISCRETE
  INT64 x;
#else
  double x;
#endif

  clear();

  s.read((char*)(&nvertex),sizeof(int));
  s.read((char*)(&background_dimension),sizeof(int));
  s.read((char*)(&euclidean),sizeof(bool));
  s.read((char*)(&relational),sizeof(bool));
  s.read((char*)(&uniform),sizeof(bool));
  s.read((char*)(&high_memory),sizeof(bool));

  if (relational) {
    for(i=0; i<nvertex; ++i) {
      for(j=1+i; j<nvertex; ++j) {
#ifdef DISCRETE
        s.read((char*)(&x),sizeof(INT64));
#else
        s.read((char*)(&x),sizeof(double));
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
        for(j=0; j<background_dimension; ++j) {
#ifdef DISCRETE
          s.read((char*)(&x),sizeof(INT64));
#else           
          s.read((char*)(&x),sizeof(double));
#endif
          xc.push_back(x);
        }
        coordinates.push_back(xc);
        xc.clear();
      }
    }
    else {
      for(i=0; i<nvertex; ++i) {
        s.read((char*)(&n),sizeof(int));
        for(j=0; j<n; ++j) {
#ifdef DISCRETE
          s.read((char*)(&x),sizeof(INT64));
#else
          s.read((char*)(&x),sizeof(double));
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
          s.read((char*)(&x),sizeof(INT64));
#else
          s.read((char*)(&x),sizeof(double));
#endif
          distances.push_back(x);
        }
      }
    }
  }
}

double Geometry::perceptual_divergence(const double* raxis,double theta,const double* translation,const double* observed) const
{
  // The method assumes that
  // a) raxis is a unit vector 
  // b) 0 <= theta < 2*pi
  // c) dimension = 3
#ifdef DEBUG
  assert(uniform && !relational && background_dimension == 3);  
#endif
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
  int i,j,k;
  std::vector<double> xc;
#ifdef DISCRETE
  std::vector<INT64> xi;
#endif
#ifdef DEBUG
  const int vsize = (signed) x.size();
#endif

  clear();

  for(i=0; i<background_dimension; ++i) {
    xc.push_back(0.0);
#ifdef DISCRETE
    xi.push_back(0);
#endif
  }

  if (unf_rnd) {
#ifdef DEBUG
    assert(2*background_dimension == vsize);
#endif
    for(i=0; i<N; ++i) {
      for(j=0; j<background_dimension; ++j) {
        xc[j] = x[2*j] + (x[2*j+1] - x[2*j])*RND.drandom();
      }
#ifdef DISCRETE
      for(j=0; j<background_dimension; ++j) {
        xi[j] = INT64(xc[j]/space_quantum);
      }
      coordinates.push_back(xi);
#else
      coordinates.push_back(xc);
#endif
    }    
  }
  else {
#ifdef DEBUG
    assert(background_dimension*N == vsize);
#endif
    for(i=0; i<N; ++i) {
      for(j=0; j<background_dimension; ++j) {
        xc[j] = x[background_dimension*i + j];
      }
#ifdef DISCRETE
      for(j=0; j<background_dimension; ++j) {
        xi[j] = INT64(xc[j]/space_quantum);
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
#ifdef _OPENMPI
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
  // Should this method operate on the assumption that we're starting from an empty geometry?
  // For the moment (December 6, 2014), we will suppose so...
  int i,j,k;
  std::vector<double> xc;
#ifdef DISCRETE
  std::vector<INT64> xi;
#endif

  clear();

  for(i=0; i<background_dimension; ++i) {
    xc.push_back(0.0);
#ifdef DISCRETE
    xi.push_back(0);
#endif
  }
  for(i=0; i<N; ++i) {
    for(j=0; j<background_dimension; ++j) {
      xc[j] = RND.nrandom(mu,sigma);
    }
#ifdef DISCRETE
    for(j=0; j<background_dimension; ++j) {
      xi[j] = INT64(xc[j]/space_quantum);
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
  // We need to set nvertex here before we start calling the compute_index
  // method

#ifdef _OPENMPI
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
  clear();

  nvertex = (signed) source.size();
#ifdef DISCRETE
  INT64 n;
  std::vector<INT64> xi;
  for(int i=0; i<nvertex; ++i) {
    for(int j=0; j<(signed) source[i].size(); ++j) {
      n = INT64(source[i][j]/space_quantum);
      xi.push_back(n);
    }
    coordinates.push_back(xi);
    xi.clear();
  }
#else
  coordinates = source;
#endif


  if (!high_memory) return;

  int i,j,k;
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
#ifdef _OPENMPI
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
  if (uniform) {
    for(int i=0; i<background_dimension; ++i) {
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
  if (relational) {

  }
  else {
    coordinates = source->coordinates;
  }
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
  if (relational) {

  }
  else {
    target->coordinates = coordinates;
  }
}

int Geometry::vertex_order(int n,int m) const
{
  if (relational || euclidean) return -1;

#ifdef DISCRETE
  INT64 sum;
#else
  double sum;
#endif

  sum = -(coordinates[n][0] - coordinates[m][0])*(coordinates[n][0] - coordinates[m][0]);
  for(int i=1; i<background_dimension; ++i) {
    sum += (coordinates[n][i] - coordinates[m][i])*(coordinates[n][i] - coordinates[m][i]);
  }
#ifdef DISCRETE
  if (sum > 0) return -1;
#else
  if (sum > 0.0) return -1;
#endif
  int output = (coordinates[n][0] < coordinates[m][0]) ? 1 : 0;
  return output;
}

void Geometry::create(int n,const std::string& type)
{
  if (!relational) {
    std::cerr << "Illegal geometric method call (initialize) for absolute model!" << std::endl;
    std::exit(1);
  }

  distances.clear();

  int i,j;
  if (type == "CARTESIAN") {
    nvertex = ipow(n,background_dimension);
    if (!high_memory) return;

    int l,p,q,s,alpha;
#ifdef DISCRETE
    INT64 delta;
#else
    double delta;
#endif
    std::vector<int> x,y;

    for(i=0; i<nvertex; ++i) {
      s = i;
      for(j=background_dimension; j>=1; j--) {
        p = ipow(n,j-1);
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
    std::cerr << "Illegal spatial type in Geometry::initialize!" << std::endl;
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
  int d1 = (signed) coordinates[n].size();
  int d2 = (signed) coordinates[m].size();
  int k=0,l = d1;
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

double Geometry::inner_product(const Matrix<double>& L,const std::vector<int>& offset) const
{
  if (relational) return 0.0;
  const int nv = L.get_nrow();
  int i,j,k,nelements;
  double value,result[nv][background_dimension],sum[background_dimension];
  double energy = 0.0;
  std::vector<double> Lx;
  std::vector<unsigned int> Lc;

  for(i=0; i<nv; ++i) {
    for(j=0; j<background_dimension; ++j) {
      sum[j] = 0.0;
    }
    L.get_row(Lx,Lc,i);
    nelements = (signed) Lx.size();
    for(j=0; j<nelements; ++j) {
      for(k=0; k<background_dimension; ++k) {
#ifdef DISCRETE
        sum[k] += space_quantum*double(coordinates[Lc[j]][k])*Lx[j];
#else
        sum[k] += coordinates[Lc[j]][k]*Lx[j];
#endif
      }
    }
    for(j=0; j<background_dimension; ++j) {
      result[i][j] = sum[j];
    }
  }

  for(i=0; i<background_dimension; ++i) {
    value = 0.0;
    for(j=0; j<nvertex; ++j) {
      if (offset[j] == -1) continue;
#ifdef DISCRETE
      value += result[offset[j]][i]*space_quantum*double(coordinates[j][i]);
#else
      value += result[offset[j]][i]*coordinates[j][i];
#endif
    }
    energy += value;
  }
  return energy;
}

void Geometry::reciprocate()
{
#ifdef DISCRETE

#else
  if (relational) {
    const int n = (signed) distances.size();
    for(int i=0; i<n; ++i) {
      distances[i] = 1.0/distances[i];
    }
  }
  else {

  }
#endif
}

void Geometry::rollback()
{
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

void Geometry::vertex_perturbation(int v)
{
  vperturb = v;
  if (relational) {
    int i,j;

    original.clear();
    for(i=0; i<nvertex; ++i) {
      if (i == v) continue;
      j = compute_index(i,v);
      original.push_back(distances[j]);
#ifdef DISCRETE
      distances[j] += INT64(RND.irandom(-50,50));
#else
      distances[j] += RND.nrandom(0.0,0.5);
#endif
    }
  }
  else {
    int k = RND.irandom(background_dimension);
    original = coordinates[v];
#ifdef DISCRETE
    coordinates[v][k] += INT64(RND.irandom(-10,10));
#else
    coordinates[v][k] += RND.nrandom(0.0,0.1);
#endif
  }
}

void Geometry::vertex_addition(const std::set<int>& antecedents)
{
  if (antecedents.empty()) {
    vertex_addition(-1);
  }
  else {
    if (RND.drandom() < 0.25) {
      vertex_addition(RND.irandom(antecedents));
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
      }
      else {
        std::vector<double> xc,avg_x;
        for(i=0; i<background_dimension; ++i) {
          avg_x.push_back(0.0);
        }
        for(it=antecedents.begin(); it!=antecedents.end(); ++it) {
          j = *it;
          for(i=0; i<background_dimension; ++i) {
#ifdef DISCRETE
            avg_x[i] += space_quantum*double(coordinates[j][i]);
#else
            avg_x[i] += coordinates[j][i];
#endif
          }
        }
        for(i=0; i<background_dimension; ++i) {
          avg_x[i] /= na;
        }
        for(i=0; i<background_dimension; ++i) {
          xc.push_back(RND.nrandom(avg_x[i],0.5));
        }
        vertex_addition(xc);
      }
    }
  }
}

void Geometry::vertex_addition(int parent,double mutation)
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
      for(i=0; i<background_dimension; ++i) {
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
      int q,p = RND.irandom(background_dimension);
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

void Geometry::mutation(int v,bool by_vertex,bool complete,double severity)
{
  int i;
#ifdef DISCRETE
  INT64 alpha = INT64((RND.nrandom(0.0,severity)/space_quantum));
#else
  double alpha = RND.nrandom(0.0,severity);
#endif

  if (by_vertex) {
    vperturb = v;
    if (relational) {
      int j;

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
      original = coordinates[v];
      if (complete) {
        for(i=0; i<(signed) coordinates[v].size(); ++i) {
#ifdef DISCRETE
          alpha = INT64((RND.nrandom(0.0,severity)/space_quantum));
#else
          alpha = RND.nrandom(0.0,severity);
#endif
          coordinates[v][i] += alpha;
        }
      }
      else {
        i = (signed) coordinates[v].size();
        coordinates[v][RND.irandom(i)] += alpha;
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
    int j,kt = 0;
    for(i=0; i<nvertex; ++i) {
      for(j=0; j<(signed) coordinates[i].size(); ++j) {
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
  if (relational || uniform) return false;
  
  int i,j,n,m;
  std::set<int> vmodified;
#ifdef DISCRETE
  std::vector<INT64> x;
#else
  std::vector<double> x;
#endif

  for(i=0; i<nvertex; ++i) {
    if (vdimension[i] == -1) continue;
    // First check to see if the dimension has changed...
    m = vdimension[i];
    n = (signed) coordinates[i].size();
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
    compute_distances(vmodified);
    return true;
  }
  return false;
}

void Geometry::compute_relational_matrices(std::vector<double>& R,std::vector<std::vector<double> >& angles) const
{
  R.clear();
  angles.clear();
  // This method supposes that we are working in Euclidean n-space with a 
  // uniform, absolute geometry
#ifdef DEBUG
  assert(background_dimension > 1 && euclidean && uniform && !relational);
#endif
  int i,j,k;
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
#ifdef _OPENMPI
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
        R[j+nvertex*i] = std::sqrt(get_distance(i,j,false));
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
#ifdef _OPENMPI
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
        r = std::sqrt(get_distance(i,j,false));
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
#ifdef _OPENMPI
#pragma omp parallel for default(shared) private(i,j,k,base,v,r,sum)
#endif
    for(i=0; i<nvertex; ++i) {
      for(j=0; j<background_dimension; ++j) {
        base[j] = coordinates[i][j];
      }
      for(j=0; j<nvertex; ++j) {
        if (i == j) continue;
        for(k=0; k<background_dimension; ++k) {
#ifdef DISCRETE
          v[k] = space_quantum*double(coordinates[j][k] - base[k]);
#else
          v[k] = coordinates[j][k] - base[k];
#endif
        }
        r = get_distance(i,j,false);
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

void Geometry::compute_distances(const std::set<int>& vmodified)
{
  if (relational || !high_memory) return;

  int i,j,k,n1,n2,in1;
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
      n1 = (signed) coordinates[i].size();
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
        n2 = (signed) coordinates[j].size();
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

void Geometry::compute_distances()
{
  if (relational || !high_memory) return;

  int i,j,k,in1;
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
#ifdef _OPENMPI
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
    int n1,n2;
#ifdef _OPENMPI
#pragma omp parallel for default(shared) private(i,j,k,n1,n2,in1,delta) schedule(dynamic,1)
#endif
    for(i=0; i<nvertex; ++i) {
      n1 = (signed) coordinates[i].size();
      in1 = i*nvertex - i*(i+1)/2;
      for(j=1+i; j<nvertex; ++j) {
        delta = pfactor*(coordinates[i][0] - coordinates[j][0])*(coordinates[i][0] - coordinates[j][0]);
        n2 = (signed) coordinates[j].size();
        n2 = (n1 <= n2) ? n1 : n2;
        for(k=1; k<n2; ++k) {
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
#ifdef DEBUG
  assert(euclidean);
#endif
  /*
  int k,in1,its = 0;
  double b,lambda,err,err_old = 0.0,delta,sum;
  std::vector<double> xnew,cdistance,BZ,Delta;
  const int M = 150;
  const double cutoff = 0.00001;
  const double pfactor = 1.0/double(nvertex);

  for(i=0; i<nvertex; ++i) {
    in1 = nvertex*i - i*(i+1)/2;
    for(j=0; j<nvertex; ++j) {
      BZ.push_back(0.0);
    }
    for(j=0; j<background_dimension; ++j) {
      x.push_back(-10.0 + 20.0*RND.drandom());
      xnew.push_back(0.0);
    }
    for(j=1+i; j<nvertex; ++j) {
      Delta[in1+j-(1+i)] = std::sqrt(distances[compute_index(i,j)]);
      cdistance.push_back(0.0);
    }
  }

  for(i=0; i<nvertex; ++i) {
    in1 = nvertex*i - i*(i+1)/2;
    for(j=1+i; j<nvertex; ++j) {
      sum = 0.0;
      for(k=0; k<background_dimension; ++k) {
        sum += (x[i*background_dimension+k]-x[j*background_dimension+k])*(x[i*background_dimension+k]-x[j*background_dimension+k]);
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
      for(j=0; j<background_dimension; ++j) {
        sum = 0.0;
        for(k=0; k<nvertex; ++k) {
          sum += BZ[i*nvertex+k]*x[background_dimension*k+j];
        }
        xnew[background_dimension*i+j] = pfactor*sum;
      }
    }
    err = 0.0;
    for(i=0; i<nvertex; ++i) {
      in1 = nvertex*i - i*(i+1)/2;
      for(j=1+i; j<nvertex; ++j) {
        sum = 0.0;
        for(k=0; k<background_dimension; ++k) {
          sum += (xnew[i*background_dimension+k]-xnew[j*background_dimension+k])*(xnew[i*background_dimension+k]-xnew[j*background_dimension+k]);
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
      D[nvertex*i+j] = get_distance(i,j,false); 
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

double SYNARMOSMA::geometry_change(const Geometry* g1,const Geometry* g2)
{
  int i,j,nva;
  bool arg1 = true;
#ifdef DISCRETE
  INT64 gdelta = 0;
#else
  double gdelta = 0.0;
#endif
#ifdef DEBUG
  assert(g1->relational == g2->relational);
  assert(g1->euclidean == g2->euclidean);
#endif
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
