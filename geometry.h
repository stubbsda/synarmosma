#include "matrix.h"

#ifndef _geometryh
#define _geometryh

namespace SYNARMOSMA {
  class Geometry;
  double geometry_change(const Geometry*,const Geometry*);

  class Geometry {
   protected:
    int nvertex;
    int vperturb;
#ifdef DISCRETE
    std::vector<INT64> original;
    std::vector<INT64> distances;
    std::vector<std::vector<INT64> > coordinates;
#else
    std::vector<double> original;
    std::vector<double> distances;
    std::vector<std::vector<double> > coordinates;
#endif
    // Whether the geometry is Euclidean or Lorentzian
    bool euclidean;
    // Whether the geometry is based on a relational or absolute model of space
    bool relational;
    // Whether the geometry is dimensionally uniform
    bool uniform;
    // Whether to fill the "distances" vector (only makes sense when relational = false)
    bool high_memory;
    // The asymptotic "flat space" dimension
    int background_dimension; 

    void clear();
    void set_default_values();
    double perceptual_divergence(const double*,double,const double*,const double*) const;
    inline int compute_index(int,int) const;
   public:
    Geometry();
    Geometry(const Geometry&);
    Geometry& operator =(const Geometry&);
    ~Geometry();
    void reciprocate();
    void initialize(bool,bool,bool,bool,int);
    int serialize(std::ofstream&) const;
    int deserialize(std::ifstream&);
    void load(const Geometry*);
    void store(Geometry*) const;
    bool consistent() const;
    void create(int,const std::string&);
    void vertex_addition(int,double mutation=1.0);
    void vertex_addition(const std::set<int>&);
    inline void vertex_addition(const std::vector<double>&);
    void multiple_vertex_addition(int,double,double);
    void multiple_vertex_addition(int,bool,const std::vector<double>&);
    void multiple_vertex_addition(const std::vector<std::vector<double> >&);    
    void get_implied_vertices(int,std::set<int>&) const;
    bool adjust_dimension(const std::vector<int>&);
    void compute_distances();
    void compute_distances(const std::set<int>&);
    double dot_product(const std::vector<double>&,const std::vector<double>&) const;
    double inner_product(const Matrix<double>&,const std::vector<int>&) const;
    int compute_coordinates(std::vector<double>&) const;
    void compute_relational_matrices(std::vector<double>&,std::vector<std::vector<double> >&) const;
    int vertex_order(int,int) const;
    void vertex_difference(int,int,std::vector<double>&) const;
    void mutation(int,bool,bool,double);
    void vertex_perturbation(int);
    void rollback();
    void geometry_restoration();
    inline bool get_euclidean() const {return euclidean;};
    inline bool get_relational() const {return relational;};
    inline bool get_uniform() const {return uniform;};
    inline bool get_memory_type() const {return high_memory;};
    inline void add(int,double);
    inline int dimension() const {return background_dimension;};
    inline Relation get_temporal_order(int,int) const;
    inline void set_element(int,double);
    inline double get_argument(const std::vector<double>&,const std::vector<double>&) const;
    inline double get_element(int) const;
    inline double get_computed_distance(int,int,bool) const;
    inline double get_distance(int,int,bool) const;
    inline double get_distance(int,const std::vector<double>&,bool) const;
    inline void get_coordinates(int,std::vector<double>&) const;
    inline void set_coordinates(int,const std::vector<double>&);
    friend double geometry_change(const Geometry*,const Geometry*);
  };

  inline int Geometry::compute_index(int v1,int v2) const
  {
#ifdef DEBUG
    assert(v1 != v2);
#endif
    int n;
    if (v1 < v2) {
      n = nvertex*v1 - v1*(1+v1)/2;
      n += (v2 - (1+v1));
    }
    else {
      n = nvertex*v2 - v2*(1+v2)/2;
      n += (v1 - (1+v2));
    }
#ifdef DEBUG
    assert(n >= 0 && n < (signed) distances.size());
#endif
    return n;
  }

  Relation Geometry::get_temporal_order(int u,int v) const 
  {
    if (relational || euclidean) return Relation::disparate;
    Relation rho = (coordinates[u][0] < coordinates[v][0]) ? Relation::before : Relation::after;
    return rho;
  }

  double Geometry::get_argument(const std::vector<double>& vx,const std::vector<double>& vy) const
  {
    double d,t,alpha,nv1 = norm(vx),nv2 = norm(vy);

    d = dot_product(vx,vy);
    alpha = d/(nv1*nv2);
    t = std::abs(alpha);
    if (uniform) {
      if (t > 1.0) alpha = alpha/t;
    }
    else {
      if (t > 1.0) alpha = 1.0/alpha;
    }
    return alpha;
  }

  double Geometry::get_element(int n) const
  {
    if (relational) return distances[n];

    int i,j,kt = 0;
    for(i=0; i<nvertex; ++i) {
      for(j=0; j<(signed) coordinates[i].size(); ++j) {
        if (kt == n) return coordinates[i][j];
        kt++;
     }
    }   
    std::cerr << "Error: Geometry element not found!" << std::endl;
    std::exit(1);
  }

  void Geometry::set_element(int n,double alpha)
  {
    if (relational) distances[n] = alpha;

    int i,j,kt = 0;
    for(i=0; i<nvertex; ++i) {
      for(j=0; j<(signed) coordinates[i].size(); ++j) {
        if (kt == n) {
          coordinates[i][j] = alpha;
          return;
        }
        kt++;
      }
    }
  }

  void Geometry::add(int n,double alpha)
  {
    if (relational) distances[n] += alpha;

    int i,j,kt = 0;
    for(i=0; i<nvertex; ++i) {
      for(j=0; j<(signed) coordinates[i].size(); ++j) {
        if (kt == n) {
          coordinates[i][j] += alpha;
          return;
        }
        kt++;
      }
    }
  }

  void Geometry::vertex_addition(const std::vector<double>& x)
  {
    if (relational) {
      std::cerr << "Illegal geometric method call for relational model!" << std::endl;
      std::exit(1);
    }
#ifdef DEBUG
    assert((signed) x.size() >= background_dimension);
#endif

#ifdef DISCRETE
    int ulimit = (uniform) ? background_dimension : (signed) x.size();
    INT64 n;
    std::vector<INT64> xt;
    for(int i=0; i<ulimit; ++i) {
      n = INT64(x[i]/space_quantum);
      xt.push_back(n);
    }
    coordinates.push_back(xt);
#else
    if (uniform) {
      std::vector<double> xt;
      for(int i=0; i<background_dimension; ++i) {
        xt.push_back(x[i]);
      }
      coordinates.push_back(xt);
    }
    else {
      coordinates.push_back(x);
    }
#endif

    if (!high_memory) {
      nvertex++;
      return;
    }

    int i,j,k = 0;
#ifdef DISCRETE
    INT64 delta;
    std::vector<INT64> ndistances;
    const int pfactor = (euclidean) ? 1 : -1;
#else
    double delta;
    std::vector<double> ndistances;
    const double pfactor = (euclidean) ? 1.0 : -1.0;
#endif

    if (uniform) {
      for(i=0; i<nvertex; ++i) {
        for(j=1+i; j<nvertex; ++j) {
          ndistances.push_back(distances[k]);
          k++;
        }
        delta = pfactor*(coordinates[i][0] - coordinates[nvertex][0])*(coordinates[i][0] - coordinates[nvertex][0]);
        for(j=1; j<background_dimension; ++j) {
          delta += (coordinates[i][j] - coordinates[nvertex][j])*(coordinates[i][j] - coordinates[nvertex][j]);
        }
        ndistances.push_back(delta);
      }
    }
    else {
      int n2,n1 = (signed) x.size();
      for(i=0; i<nvertex; ++i) {
        for(j=1+i; j<nvertex; ++j) {
          ndistances.push_back(distances[k]);
          k++;
        }
        delta = pfactor*(coordinates[i][0] - coordinates[nvertex][0])*(coordinates[i][0] - coordinates[nvertex][0]);
        n2 = (signed) coordinates[i].size();
        n2 = (n2 > n1) ? n1 : n2;
        for(j=1; j<n2; ++j) {
          delta += (coordinates[i][j] - coordinates[nvertex][j])*(coordinates[i][j] - coordinates[nvertex][j]);
        }
        ndistances.push_back(delta);
      }
    }
    distances = ndistances;
    nvertex++;
  }

  void Geometry::get_coordinates(int v,std::vector<double>& x) const
  {
    if (relational) {
      std::cerr << "Illegal geometric method call for relational model!" << std::endl;
      std::exit(1);
    }
#ifdef DISCRETE
    x.clear();
    for(int i=0; i<(signed) coordinates[v].size(); ++i) {
      x.push_back(space_quantum*double(coordinates[v][i]));
    }
#else
    x = coordinates[v];
#endif
  }

  void Geometry::set_coordinates(int v,const std::vector<double>& x)
  {
    if (relational) {
      std::cerr << "Illegal geometric method call for relational model!" << std::endl;
      std::exit(1);
    }
#ifdef DEBUG
    assert((signed) x.size() >= background_dimension);
#endif

#ifdef DISCRETE
    int ulimit = (uniform) ? background_dimension : (signed) x.size();
    INT64 n;
    std::vector<INT64> xt;
    for(int i=0; i<ulimit; ++i) {
      n = INT64(x[i]/space_quantum);
      xt.push_back(n);
    }
    coordinates[v] = xt;
#else
    if (uniform) {
      std::vector<double> xt;
      for(int i=0; i<background_dimension; ++i) {
        xt.push_back(x[i]);
      }
      coordinates[v] = xt;
    }
    else {
      coordinates[v] = x;
    }
#endif
  }

  double Geometry::get_distance(int v,const std::vector<double>& x,bool lorentzian) const 
  {
    if (relational) {
      std::cerr << "Illegal geometric method call for relational model!" << std::endl;
      std::exit(1);
    }
#ifdef DEBUG
    assert((signed) x.size() >= background_dimension);
#endif

#ifdef DISCRETE    
    std::vector<double> base;
    for(int i=0; i<(signed) coordinates[v].size(); ++i) {
      base.push_back(space_quantum*double(coordinates[v][i]));
    }
#else
    std::vector<double> base = coordinates[v];
#endif
    double delta = (base[0] - x[0])*(base[0] - x[0]);
    if (lorentzian) delta = -delta;

    if (uniform) {
      for(int i=1; i<background_dimension; ++i) {
        delta += (base[i] - x[i])*(base[i] - x[i]);
      }
    }
    else {
      int n1 = (signed) base.size();
      int n2 = (signed) x.size();
      n1 = (n1 <= n2) ? n1 : n2;
      for(int i=1; i<n1; ++i) {
        delta += (base[i] - x[i])*(base[i] - x[i]);
      }
    }
    return delta;
  }

  double Geometry::get_distance(int v1,int v2,bool lorentzian) const
  {
    if (!high_memory) return get_computed_distance(v1,v2,lorentzian);
#ifdef DEBUG    
    assert(v1 != v2);
#endif
    int n = compute_index(v1,v2);  

#ifdef DISCRETE
    double l = space_quantum*space_quantum*double(distances[n]);
#else
    double l = distances[n];
#endif
    if (!lorentzian && l < 0.0) l = -l;
    return l;
  }

  double Geometry::get_computed_distance(int v1,int v2,bool lorentzian) const
  {
#ifdef DEBUG
    assert(v1 != v2);
#endif
    double l;
    int n;
    if (relational) {
#ifdef DISCRETE
      l = space_quantum*double(distances[compute_index(v1,v2)]);
#else
      l = distances[compute_index(v1,v2)];
#endif
      if (!lorentzian && l < 0.0) l = -l;
    }
    else {
      int m;

      if (uniform) {
        m = background_dimension;
      }
      else {
        m = (signed) coordinates[v1].size();
        n = (signed) coordinates[v2].size();
        m = (m <= n) ? m : n;
      }
#ifdef DISCRETE
      INT64 q = (coordinates[v1][0] - coordinates[v2][0])*(coordinates[v1][0] - coordinates[v2][0]);
      if (lorentzian) q = -q;
      for(n=1; n<m; ++n) {
        q += (coordinates[v1][n] - coordinates[v2][n])*(coordinates[v1][n] - coordinates[v2][n]);
      }
      l = space_quantum*space_quantum*double(q);
#else
      l = (coordinates[v1][0] - coordinates[v2][0])*(coordinates[v1][0] - coordinates[v2][0]);
      if (lorentzian) l = -l;
      for(n=1; n<m; ++n) {
        l += (coordinates[v1][n] - coordinates[v2][n])*(coordinates[v1][n] - coordinates[v2][n]);
      }
#endif
    }
    return l;
  }
}
#endif


