#include "random.h"
#include "matrix.h"

#ifndef _geometryh
#define _geometryh

namespace SYNARMOSMA {
  class Geometry;
  double geometry_change(const Geometry*,const Geometry*);

  /// A class representing the geometry of a collection of points, stored either as coordinates or in relational form.
  class Geometry {
   protected:
    /// Ths number of points or vertices whose geometry is 
    /// being described by this instance of the class.
    int nvertex = 0;
    int vperturb = -1;
#ifdef DISCRETE
    std::vector<INT64> original;
    std::vector<INT64> distances;
    std::vector<std::vector<INT64> > coordinates;
#else
    std::vector<double> original;
    /// This vector property contains all of the (nvertex-1)*nvertex/2 distances between the 
    /// points and is used when Geometry::high_memory is true to speed up geometry calculations 
    /// by precomputing all the distances. 
    std::vector<double> distances;
    /// This property stores the coordinates of each point, as a distinct vector, and of course 
    /// is meaningful only when Geometry::relational is false. The number of coordinates for each 
    /// point may not be the same, unless Geometry::uniform is true.
    std::vector<std::vector<double> > coordinates;
#endif
    /// This Boolean property determines whether the geometry is 
    /// Euclidean or Lorentzian.
    bool euclidean = true;
    /// This Boolean property determines whether the geometry is 
    /// based on a relational or absolute (coordinate-based) model 
    /// of space.
    bool relational = false;
    /// This Boolean property determines whether the geometry is 
    /// dimensionally uniform, i.e. each point has the same dimensionality.
    bool uniform = true;
    /// This Boolean property determines whether to fill the Geometry::distances 
    /// vector and is only meaningful when Geometry::relational is false.
    bool high_memory = true;
    /// This non-negative property corresponds to the asymptotic "flat 
    /// space" dimension of the space.
    unsigned int background_dimension = 3; 

    /// This method clears the vector properties and sets the scalar class properties to their default value.
    void clear();
    double perceptual_divergence(const double*,double,const double*,const double*) const;
    inline int compute_index(int,int) const;
   public:
    /// The default constructor which does nothing.
    Geometry();
    /// The copy constructor that copies over all the class properties.
    Geometry(const Geometry&);
    /// The overloaded assignment operator that sets all the class properties to those of the source.
    Geometry& operator =(const Geometry&);
    /// The destructor which does nothing.
    ~Geometry();
    void initialize(bool,bool,bool,bool,int);
    /// This method writes the instance properties to a binary disk file and returns the number of bytes written to the file.
    int serialize(std::ofstream&) const;
    /// This method calls the clear() method on the instance and then reads the properties from a binary disk file and returns the number of bytes read.
    int deserialize(std::ifstream&);
    void load(const Geometry*);
    void store(Geometry*) const;
    bool consistent() const;
    void create(int,std::string&);
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
    unsigned int compute_coordinates(std::vector<double>&) const;
    void compute_relational_matrices(std::vector<double>&,std::vector<std::vector<double> >&) const;
    int vertex_order(int,int) const;
    void vertex_difference(int,int,std::vector<double>&) const;
    void mutation(int,bool,bool,double);
    void vertex_perturbation(int);
    void rollback();
    void geometry_restoration();
    /// This method returns the value of the Geometry::euclidean property.
    inline bool get_euclidean() const {return euclidean;};
    /// This method returns the value of the Geometry::relational property.
    inline bool get_relational() const {return relational;};
    /// This method returns the value of the Geometry::uniform property.
    inline bool get_uniform() const {return uniform;};
    /// This method returns the value of the Geometry::high_memory property.
    inline bool get_memory_type() const {return high_memory;};
    inline void add(int,double);
    inline unsigned int dimension() const {return background_dimension;};
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
    if (v1 == v2) throw std::invalid_argument("The vertex arguments in Geometry::compute_index must be distinct!");
    if (v1 < 0 || v1 >= nvertex) throw std::invalid_argument("The first vertex argument in Geometry::compute_index has an illegal value!");
    if (v2 < 0 || v2 >= nvertex) throw std::invalid_argument("The second vertex argument in Geometry::compute_index has an illegal value!");

    int n;
    if (v1 < v2) {
      n = nvertex*v1 - v1*(1+v1)/2;
      n += (v2 - (1+v1));
    }
    else {
      n = nvertex*v2 - v2*(1+v2)/2;
      n += (v1 - (1+v2));
    }
    if (n < 0 || n >= (signed) distances.size()) throw std::runtime_error("The computed index value in Geometry::compute_index is illegal!");

    return n;
  }

  Relation Geometry::get_temporal_order(int u,int v) const 
  {
    auto rho = Relation::disparate;
    if (relational || euclidean) return rho;
    // Check that the interval is timelike...
    if (get_distance(u,v,true) < -std::numeric_limits<double>::epsilon()) {
      rho = (coordinates[u][0] < coordinates[v][0]) ? Relation::before : Relation::after;
    }
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
    throw std::invalid_argument("Missing element in Geometry::get_element method!");
  }

  void Geometry::set_element(int n,double alpha)
  {
    if (relational) distances[n] = alpha;

    int i,kt = 0;
    unsigned int j;
    for(i=0; i<nvertex; ++i) {
      for(j=0; j<coordinates[i].size(); ++j) {
        if (kt == n) {
          coordinates[i][j] = alpha;
          return;
        }
        kt++;
      }
    }
    throw std::invalid_argument("Missing element in Geometry::set_element method!");
  }

  void Geometry::add(int n,double alpha)
  {
    if (relational) distances[n] += alpha;

    int i,kt = 0;
    unsigned int j;
    for(i=0; i<nvertex; ++i) {
      for(j=0; j<coordinates[i].size(); ++j) {
        if (kt == n) {
          coordinates[i][j] += alpha;
          return;
        }
        kt++;
      }
    }
    throw std::invalid_argument("Missing element in Geometry::add method!");
  }

  void Geometry::vertex_addition(const std::vector<double>& x)
  {
    if (relational) throw std::runtime_error("Illegal method call (Geometry::vertex_addition) for relational model!");
    if (x.size() < background_dimension) throw std::invalid_argument("The length of the vector argument of Geometry::vertex_addition must not be less than the background dimension!");

    unsigned int i;
#ifdef DISCRETE
    unsigned int ulimit = (uniform) ? background_dimension : x.size();
    INT64 n;
    std::vector<INT64> xt;
    for(i=0; i<ulimit; ++i) {
      n = INT64(x[i]/space_quantum);
      xt.push_back(n);
    }
    coordinates.push_back(xt);
#else
    if (uniform) {
      std::vector<double> xt;
      for(i=0; i<background_dimension; ++i) {
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

    unsigned int j,k = 0;
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
      for(int l=0; l<nvertex; ++l) {
        for(int p=1+l; p<nvertex; ++p) {
          ndistances.push_back(distances[k]);
          k++;
        }
        delta = pfactor*(coordinates[l][0] - coordinates[nvertex][0])*(coordinates[l][0] - coordinates[nvertex][0]);
        for(j=1; j<background_dimension; ++j) {
          delta += (coordinates[l][j] - coordinates[nvertex][j])*(coordinates[l][j] - coordinates[nvertex][j]);
        }
        ndistances.push_back(delta);
      }
    }
    else {
      unsigned int n2,n1 = x.size();
      for(int l=0; l<nvertex; ++l) {
        for(int p=1+l; p<nvertex; ++p) {
          ndistances.push_back(distances[k]);
          k++;
        }
        delta = pfactor*(coordinates[l][0] - coordinates[nvertex][0])*(coordinates[l][0] - coordinates[nvertex][0]);
        n2 = coordinates[i].size();
        n2 = (n2 > n1) ? n1 : n2;
        for(j=1; j<n2; ++j) {
          delta += (coordinates[l][j] - coordinates[nvertex][j])*(coordinates[l][j] - coordinates[nvertex][j]);
        }
        ndistances.push_back(delta);
      }
    }
    distances = ndistances;
    nvertex++;
  }

  void Geometry::get_coordinates(int v,std::vector<double>& x) const
  {
    if (relational) throw std::runtime_error("Illegal method call (Geometry::get_coordinates) for relational model!");
#ifdef DISCRETE
    x.clear();
    unsigned int i;
    for(i=0; i<coordinates[v].size(); ++i) {
      x.push_back(space_quantum*double(coordinates[v][i]));
    }
#else
    x = coordinates[v];
#endif
  }

  void Geometry::set_coordinates(int v,const std::vector<double>& x)
  {
    if (relational) throw std::runtime_error("Illegal method call (Geometry::set_coordinates) for relational model!");
    if (x.size() < background_dimension) throw std::invalid_argument("The length of the vector argument of Geometry::set_coordinates must not be less than the background dimension!");

    unsigned int i;
#ifdef DISCRETE
    unsigned int ulimit = (uniform) ? background_dimension : x.size();
    INT64 n;
    std::vector<INT64> xt;
    for(i=0; i<ulimit; ++i) {
      n = INT64(x[i]/space_quantum);
      xt.push_back(n);
    }
    coordinates[v] = xt;
#else
    if (uniform) {
      std::vector<double> xt;
      for(i=0; i<background_dimension; ++i) {
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
    if (relational) throw std::runtime_error("Illegal method call (Geometry::get_distance) for relational model!");
    if (x.size() < background_dimension) throw std::invalid_argument("The length of the vector argument of Geometry::get_distance must not be less than the background dimension!");

    unsigned int i;
#ifdef DISCRETE    
    std::vector<double> base;
    for(i=0; i<coordinates[v].size(); ++i) {
      base.push_back(space_quantum*double(coordinates[v][i]));
    }
#else
    std::vector<double> base = coordinates[v];
#endif
    double delta = (base[0] - x[0])*(base[0] - x[0]);
    if (lorentzian) delta = -delta;

    if (uniform) {
      for(i=1; i<background_dimension; ++i) {
        delta += (base[i] - x[i])*(base[i] - x[i]);
      }
    }
    else {
      unsigned int n1 = base.size(),n2 = x.size();

      n1 = (n1 <= n2) ? n1 : n2;
      for(i=1; i<n1; ++i) {
        delta += (base[i] - x[i])*(base[i] - x[i]);
      }
    }
    return delta;
  }

  double Geometry::get_distance(int v1,int v2,bool lorentzian) const
  {
    if (v1 == v2) throw std::invalid_argument("The vertex arguments in Geometry::get_distance must be distinct!");

    if (!high_memory) return get_computed_distance(v1,v2,lorentzian);

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
    if (v1 == v2) throw std::invalid_argument("The vertex arguments in Geometry::get_computed_distance must be distinct!");

    double l;
    unsigned int n,m;
    if (relational) {
#ifdef DISCRETE
      l = space_quantum*double(distances[compute_index(v1,v2)]);
#else
      l = distances[compute_index(v1,v2)];
#endif
      if (!lorentzian && l < 0.0) l = -l;
    }
    else {
      if (uniform) {
        m = background_dimension;
      }
      else {
        m = coordinates[v1].size();
        n = coordinates[v2].size();
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


