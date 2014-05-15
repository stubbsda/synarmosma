#include "global.h"

class Geometry {
 private:
  int nvertex;
  int vperturb;
  std::vector<double> original;
  std::vector<double> distances;
#ifdef LEIBNIZ
  hash_map index;
#else
  std::vector<std::vector<double> > coordinates;
#endif

 public:
  // Whether the geometry is Euclidean or Lorentzian
  bool euclidean;
  // The asymptotic "flat space" dimension
  static const int background_dimension = 2;

  Geometry();
  Geometry(bool);
  Geometry(const Geometry&);
  Geometry& operator =(const Geometry&);
  ~Geometry();
  void reciprocate();
  void clear();
  void serialize(std::ofstream&) const;
  void deserialize(std::ifstream&);
  void load(const Geometry*);
  void store(Geometry*) const;
  void initialize(int,const std::string&);
  void add_vertex(int,double mutation=1.0);
  void add_vertex(const std::set<int>&);
  inline void add_vertex(const std::vector<double>&);
  void get_implied_vertices(int,std::set<int>&) const;
  bool adjust_dimension(const std::vector<int>&);
  void compute_distances();
  void compute_distances(const std::set<int>&);
  double dot_product(const std::vector<double>&,const std::vector<double>&) const;
  double inner_product(const std::vector<double>*,const std::vector<int>&,int) const;
  int compute_coordinates(std::vector<double>&) const;
  int vertex_order(int,int) const;
  void vertex_difference(int,int,std::vector<double>&) const;
  void additive_modification(int,bool,double,double);
  void multiplicative_modification(int,bool,double,double);
  void perturb_vertex(int);
  void rollback();
  void geometry_modification(int,double,double);
  void geometry_restoration();
  void set_metric(const std::string&);
  inline void set_metric(bool t) {euclidean = t;};
  inline void add(int,double);
  inline void set_element(int,double);
  inline double get_element(int) const;
  inline double get_computed_distance(int,int,bool) const;
  inline double get_distance(int,int,bool) const;
  inline double get_distance(int,const std::vector<double>&,bool) const;
  inline void get_coordinates(int,std::vector<double>&) const;
  inline void set_coordinates(int,const std::vector<double>&);
  friend double geometry_change(const Geometry*,const Geometry*);
};

double Geometry::get_element(int n) const
{
#ifdef LEIBNIZ
  return distances[n];
#else
  int i,j,kt = 0;
  for(i=0; i<nvertex; ++i) {
    assert(coordinates[i].size() == 2);
    for(j=0; j<(signed) coordinates[i].size(); ++j) {
      if (kt == n) return coordinates[i][j];
      kt++;
    }
  }
  std::cerr << "Error: Geometry element not found!" << std::endl;
  std::exit(1);
#endif
}

void Geometry::set_element(int n,double alpha)
{
#ifdef LEIBNIZ
  distances[n] = alpha;
#else
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
#endif
}

void Geometry::add(int n,double alpha)
{
#ifdef LEIBNIZ
  distances[n] += alpha;
#else
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
#endif
}

void Geometry::add_vertex(const std::vector<double>& x)
{
#ifdef LEIBNIZ
  std::cerr << "Illegal geometric method call for relational model!" << std::endl;
  std::exit(1);
#else
  coordinates.push_back(x);
  nvertex++;
#endif
}

void Geometry::get_coordinates(int v,std::vector<double>& x) const
{
#ifdef LEIBNIZ
  std::cerr << "Illegal geometric method call for relational model!" << std::endl;
  std::exit(1);
#else
  x = coordinates[v];
#endif
}

void Geometry::set_coordinates(int v,const std::vector<double>& x)
{
#ifdef LEIBNIZ
  std::cerr << "Illegal geometric method call for relational model!" << std::endl;
  std::exit(1);
#else
  coordinates[v] = x;
#endif
}

double Geometry::get_distance(int v,const std::vector<double>& x,bool causal) const 
{
  double delta = 0.0;
#ifdef LEIBNIZ
  std::cerr << "Illegal geometric method call for relational model!" << std::endl;
  std::exit(1);
#else
  if (causal) {
    delta = -(coordinates[v][0] - x[0])*(coordinates[v][0] - x[0]);
    for(int i=1; i<background_dimension; ++i) {
      delta += (coordinates[v][i] - x[i])*(coordinates[v][i] - x[i]);
    }
  }
  else {
    for(int i=0; i<background_dimension; ++i) {
      delta += (coordinates[v][i] - x[i])*(coordinates[v][i] - x[i]);
    }
  }
#endif
  return delta;
}

double Geometry::get_distance(int v1,int v2,bool lorentzian) const
{
  assert(v1 != v2);
  int n;
  double l;
#ifdef LEIBNIZ
  hash_map::const_iterator qt = index.find(make_key(v1,v2));
  n = qt->second;
#else
  if (v1 < v2) {
    n = nvertex*v1 - v1*(1+v1)/2;
    n += (v2 - (1+v1));
  }
  else {
    n = nvertex*v2 - v2*(1+v2)/2;
    n += (v1 - (1+v2));
  }
#endif
  l = distances[n];
  if (!lorentzian && l < 0.0) l = -l;
  return l;
}

double Geometry::get_computed_distance(int v1,int v2,bool lorentzian) const
{
  assert(v1 != v2);
  double l;
  int n;
#ifdef LEIBNIZ
  hash_map::const_iterator qt = index.find(make_key(v1,v2));
  assert(qt != index.end());
  n = qt->second;
  l = distances[n];
  if (!lorentzian && l < 0.0) l = -l;
#else
  l = (coordinates[v1][0] - coordinates[v2][0])*(coordinates[v1][0] - coordinates[v2][0]);
  if (lorentzian) l = -l;
  for(n=1; n<background_dimension; ++n) {
    l += (coordinates[v1][n] - coordinates[v2][n])*(coordinates[v1][n] - coordinates[v2][n]);
  }
#endif
  return l;
}


