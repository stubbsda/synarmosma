#include "random.h"
#include "matrix.h"

#ifndef _geometryh
#define _geometryh

namespace SYNARMOSMA {
  template<class kind>
  class Geometry;

  template<class kind>
  double geometry_change(const Geometry<kind>*,const Geometry<kind>*);

  /// A class representing the geometry of a collection of points, stored either as coordinates or in relational form.
  template<class kind>
  class Geometry {
   protected:
    /// This property represents the number of points or vertices whose geometry is 
    /// being described by this instance of the class.
    int nvertex = 0;
    /// This property is used to store the index of the point whose geometry has been 
    /// perturbed and allows the geometry to be changed in an atomic manner, with the 
    /// possibility of undoing such a change and without having to recalculate all of 
    /// inter-vertex distances. 
    int vperturb = -1;
    /// This vector property is like Geometry::vperturb used to store the distance values 
    /// of the vertex whose geometry is perturbed so that they can be restored in case the 
    /// change is rolled back.
    std::vector<kind> original;
    /// This vector property contains all of the (nvertex-1)*nvertex/2 squared distances between 
    /// the points. It is used when Geometry::relational is true or when Geometry::high_memory is 
    /// true to speed up geometry calculations by precomputing all the squared distances. 
    std::vector<kind> distances;
    /// This property stores the coordinates of each point, as a distinct vector, and of course 
    /// is meaningful only when Geometry::relational is false. The number of coordinates for each 
    /// point may not be the same, unless Geometry::uniform is true.
    std::vector<std::vector<kind> > coordinates;
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
    /// This constant represents the smallest possible spatial separation between two vertices.
    static const double space_quantum; 

    /// This method accepts as its input an axis of rotation (the first argument), an angle (second argument), a translation vector (third argument) and finally a set of observational locations (the final argument), the perceived three-dimensional coordinates for each vertex. The method applies the rotation and translation to the coordinates of every vertex in the geometry and then computes the distance from this value to the observed vertex location, adds together these distances and returns their arithmetic mean.
    double perceptual_divergence(const double*,double,const double*,const double*) const;
    /// This method calculates the index into the vector Geometry::distances given the indices of two vertices, the method's arguments. 
    int compute_index(int,int) const;
    /// This method converts the first argument to a double and returns this, by multiplying it by Geometry::space_quantum raised to the power of the second argument.
    double convert(kind,int = 1) const;
    /// This method converts its argument to the base type of this class and returns this, by first dividing it by Geometry::space_quantum and then truncating the fractional part.
    kind invert(double) const;
   public:
    /// The default constructor which does nothing.
    Geometry();
    /// The copy constructor that copies over all the class properties.
    Geometry(const Geometry&);
    /// The overloaded assignment operator that sets all the class properties to those of the source.
    Geometry& operator =(const Geometry&);
    /// The destructor which does nothing.
    ~Geometry();
    /// This method calls clear() and then initializes the scalar properties of this class: Geometry::euclidean, Geometry::relational, Geometry::uniform, Geometry::high_memory and Geometry::background_dimension.
    void initialize(bool,bool,bool,bool,int);
    /// This method writes the instance properties to a binary disk file and returns the number of bytes written to the file.
    int serialize(std::ofstream&) const;
    /// This method calls the clear() method on the instance and then reads the properties from a binary disk file and returns the number of bytes read.
    int deserialize(std::ifstream&);
    /// This method clears the vector properties and sets the scalar class properties to their default value.
    void clear();
    /// This method loads the properties of the method's argument to the current instance.
    void load(const Geometry*);
    /// This method stores the properties of the current instance in the method's argument.
    void store(Geometry*) const;
    /// This method verifies that the lengths of Geometry::coordinates and Geometry::distances and then tests these two vectors for any NaN elements; if there are any problems it returns false.
    bool consistent() const;
    /// This method creates a relational geometry involving n vertices (where n is the first argument) and a geometric model, described in the second argument: CARTESIAN, SINGLETON, MONOPLEX or RANDOM; by default we suppose a Cartesian geometric model.
    void create(int,const std::string& = std::string("CARTESIAN"));
    /// This method adds a vertex to the geometry; the first argument is the vertex's parent and the second its distance from this parent. If the first argument is -1, the vertex is placed randomly and the second argument is ignored. The return value is the index of the new vertex.
    int vertex_addition(int,double mutation=1.0);
    /// This method adds a vertex to the geometry; its unique argument is a set of parents for this new vertex. If the set is empty the vertex is placed randomly, otherwise its localization is a random Gaussian variate whose mean is the average of the location of its parents. The return value is the index of the new vertex.
    int vertex_addition(const std::set<int>&);
    /// This method adds a vertex to the geometry; its unique argument is the coordinates for the new vertex and it should only be called therefore when Geometry::relational is false. The return value is the index of the new vertex.
    int vertex_addition(const std::vector<double>&);
    /// This method adds multiple vertices to the geometry - it begins by calling the clear() method so that it starts from an empty geometry, then adds n vertices where n is the first argument. The second and third arguments are the mean and standard deviation for the coordinates of these n vertices; note that this method is not compatible with a relational geometry.
    void multiple_vertex_addition(int,double,double);
    /// This method adds multiple vertices to the geometry - it begins by calling the clear() method so that it starts from an empty geometry, then adds n vertices where n is the first argument. The second argument determines whether or not the coordinates are drawn from a uniform random distribution - if true the third argument must be a vector of length twice the Geometry::background_dimension containing the coordinate limits, if false then the third argument must contain all of the coordinates (and so be of length n times Geometry::background_dimension). Note that this method is not compatible with a relational geometry.
    void multiple_vertex_addition(int,bool,const std::vector<double>&);
    /// This method adds multiple vertices to the geometry - it begins by calling the clear() method so that it starts from an empty geometry, then adds n vertices where n is the length of the outer vector of its unique argument. This argument is the set of coordinates for the vertices in this geometry; note that this method is not compatible with a relational geometry.
    void multiple_vertex_addition(const std::vector<std::vector<double> >&);
    /// This method determines which vertices correspond to the index that is the method's first argument and write these vertices in the second argument. If Geometry::relational is true these are the two vertices whose distance is indexed, otherwise it is the single vertex whose coordinate index is indicated.    
    void get_implied_vertices(int,std::set<int>&) const;
    /// This method alters the size of the tuple associated with a vertex depending on its dimension, which is contained in the method's argument. If any changes have to be made to the geometry then the method returns true and false otherwise.
    bool adjust_dimension(const std::vector<int>&);
    /// This method computes the squared distances based on the coordinate values; if Geometry::relational is true or Geometry::high_memory is false, the method returns immediately.
    void compute_squared_distances();
    /// This method computes the squared distances between the vertices in the method's argument and the rest of the vertices in the geometry; if Geometry::relational is true or Geometry::high_memory is false, the method returns immediately.
    void compute_squared_distances(const std::set<int>&);
    /// This method computes and returns the dot product of the two vectors that are its argument, according to the signature - Euclidean or Lorentzian - of the geometry; if the vectors do not have the same length, the method pads the shorter vector with elements that random Gaussian variates with mean zero and a shrinking standard deviation.
    double dot_product(const std::vector<double>&,const std::vector<double>&) const;
    /// This method calculates the coordinates of the geometry and writes them to the method's argument. When Geometry::relational is false this method has relatively little to do but when it is true, it uses multi-dimensional scaling to compute coordinates for the vertices. The method returns the maximum dimension of the vertices.
    unsigned int compute_coordinates(std::vector<double>&) const;
    /// THis method calculates a relational representation of the coordinate geometry, contained in its two arguments - the first is the vector of distances among the vertices and the second is a matrix of angles among them. The method is only meaningful when Geometry::relational is false, Geometry::background_dimension is greater than one and both Geometry::uniform and Geometry::euclidean are true.
    void compute_relational_matrices(std::vector<double>&,std::vector<std::vector<double> >&) const;
    /// This method calculates the difference between the coordinate vectors for the two vertices that are its first arguments and writes the result in the final argument. The method is obviously meaningless in a relational geometry and pads the coordinate vectors as necessary with random Gaussian variates if the geometry isn't uniform.
    void vertex_difference(int,int,std::vector<double>&) const;
    /// This method perturbs the geometry of the vertex whose index is the method's first argument. The second and third arguments determine if the mutation is by vertex and if it is complete (i.e. affecting all coordinate dimensions) while the final argument is the severity. The mutation is represented as the addition of a random Gaussian variate with mean zero and a small standard deviation to the existing value, while backing up the former value using the Geometry::vperturb and Geometry::original properties.
    void mutation(int,bool,bool,double);
    /// If this method is called, then the geometry is restored to its original nature, based on the value of Geometry::vperturb and Geometry::original, after a call to vertex_perturbation(). The argument determines whether or not the rollback is minimal, based on the severity of the perturbation of the geometry.
    void rollback(bool);
    /// This method returns the value of the Geometry::euclidean property.
    bool get_euclidean() const;
    /// This method returns the value of the Geometry::relational property.
    bool get_relational() const;
    /// This method returns the value of the Geometry::uniform property.
    bool get_uniform() const;
    /// This method returns the value of the Geometry::high_memory property.
    bool get_memory_type() const;
    /// This method returns the value the Geometry::background_dimension property.
    unsigned int dimension() const;
    /// This method adds the second argument to the value of Geometry::distances[n] (if the geometry is relational) or to the n-th element in Geometry::coordinates (if not).
    void add(int,double);
    /// This method returns the nature of the temporal ordering between the two points indicated by the method's arguments: DISPARATE, BEFORE or AFTER. 
    Relation get_temporal_order(int,int) const;
    /// This method computes and returns the quotient of the dot product of its two arguments, divided by the product of their norm, which in a Euclidean context is equal to the cosine of the angle between them. 
    double get_argument(const std::vector<double>&,const std::vector<double>&) const;
    /// This method gets the value of the n-th element in Geometry::coordinates by counting through each tuple from Geometry::coordinates[0][0]; if the geometry is relational it gets the value of Geometry::distances[n].
    double get_element(int) const;
    /// This method sets the value of the n-th element in Geometry::coordinates to its second argument by counting through each tuple from Geometry::coordinates[0][0]; if the geometry is relational it sets the value of Geometry::distances[n] to the second argument.
    void set_element(int,double);
    /// This method returns the squared distance between the two vertices specified by its first two arguments; if the third argument is true, then the method will recompute the squared distance (this is meaningless if Geometry::relational is true). 
    double get_squared_distance(int,int,bool) const;
    /// This method computes and returns the squared distance between the point whose index is the first argument and the coordinate vector that is the second argument; this method is meaningless if Geometry::relational is true.
    double get_squared_distance(int,const std::vector<double>&) const;
    /// This method gets the coordinate vector for the point whose index is the first argument and sets the second argument equal to this coordinate tuple.
    void get_coordinates(int,std::vector<double>&) const;
    /// This method sets the coordinate vector for the point whose index is the first argument to the method's second argument. 
    void set_coordinates(int,const std::vector<double>&);
    /// This method computes the difference between two instances of this class, returning a non-negative floating point number; the greater the number, the more substantial the difference between the two geometries.
    friend double geometry_change<>(const Geometry<kind>*,const Geometry<kind>*);
  };

  template<>
  /// This method has been specialized for the double base type, since in this case the method does nothing - it just returns its first argument.
  inline double Geometry<double>::convert(double x,int n) const
  {
    return x;
  }

  template<>
  /// This method has been specialized for the double base type, since in this case the method does nothing - it just returns its first argument.
  inline double Geometry<double>::invert(double x) const
  {
    return x;
  }

  template<class kind>
  inline double Geometry<kind>::convert(kind x,int n) const
  {
    double pfactor = space_quantum;
    for(int i=1; i<n; ++i) {
      pfactor *= space_quantum;
    }
    return pfactor*double(x);
  }

  template<class kind>
  inline kind Geometry<kind>::invert(double x) const
  {
    kind n = kind(x/space_quantum);
    return n;
  }

  template<class kind>
  inline bool Geometry<kind>::get_euclidean() const
  {
    return euclidean;
  }

  template<class kind>
  inline bool Geometry<kind>::get_relational() const
  {
    return relational;
  }

  template<class kind>
  inline bool Geometry<kind>::get_uniform() const
  {
    return uniform;
  }

  template<class kind>
  inline bool Geometry<kind>::get_memory_type() const
  {
    return high_memory;
  }

  template<class kind>
  inline unsigned int Geometry<kind>::dimension() const
  {
    return background_dimension;
  }

  template<class kind>
  inline int Geometry<kind>::compute_index(int v1,int v2) const
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

  template<class kind>
  inline Relation Geometry<kind>::get_temporal_order(int u,int v) const 
  {
    if (u == v) throw std::invalid_argument("The vertex arguments in Geometry::get_temporal_order must be distinct!");
    Relation rho = Relation::disparate;
    if (relational || euclidean) return rho;
    // Check that the interval is timelike...
    if (get_squared_distance(u,v,false) < -std::numeric_limits<double>::epsilon()) {
      rho = (coordinates[u][0] < coordinates[v][0]) ? Relation::before : Relation::after;
    }
    return rho;
  }

  template<class kind>
  inline double Geometry<kind>::get_argument(const std::vector<double>& vx,const std::vector<double>& vy) const
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

  template<class kind>
  inline double Geometry<kind>::get_element(int n) const
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

  template<class kind>
  inline void Geometry<kind>::set_element(int n,double alpha)
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

  template<class kind>
  inline void Geometry<kind>::add(int n,double alpha)
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

  template<class kind>
  inline int Geometry<kind>::vertex_addition(const std::vector<double>& x)
  {
    if (relational) throw std::runtime_error("Illegal method call (Geometry::vertex_addition) for relational model!");
    if (x.size() < background_dimension) throw std::invalid_argument("The length of the vector argument of Geometry::vertex_addition must not be less than the background dimension!");

    unsigned int i;
    std::vector<kind> xt;
    unsigned int ulimit = (uniform) ? background_dimension : x.size();
    for(i=0; i<ulimit; ++i) {
      xt.push_back(invert(x[i]));
    }
    coordinates.push_back(xt);

    if (!high_memory) {
      nvertex++;
      return nvertex;
    }

    unsigned int j,k = 0;
    kind delta;
    std::vector<kind> ndistances;
    const kind pfactor = (euclidean) ? kind(1) : kind(-1);

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
        n2 = coordinates[l].size();
        n2 = (n2 > n1) ? n1 : n2;
        for(j=1; j<n2; ++j) {
          delta += (coordinates[l][j] - coordinates[nvertex][j])*(coordinates[l][j] - coordinates[nvertex][j]);
        }
        ndistances.push_back(delta);
      }
    }
    distances = ndistances;
    nvertex++;
    return nvertex;
  }

  template<class kind>
  inline void Geometry<kind>::get_coordinates(int v,std::vector<double>& x) const
  {
    if (relational) throw std::runtime_error("Illegal method call (Geometry::get_coordinates) for relational model!");
    unsigned int i;

    x.clear();

    for(i=0; i<coordinates[v].size(); ++i) {
      x.push_back(convert(coordinates[v][i]));
    }
  }

  template<class kind>
  inline void Geometry<kind>::set_coordinates(int v,const std::vector<double>& x)
  {
    if (relational) throw std::runtime_error("Illegal method call (Geometry::set_coordinates) for relational model!");
    if (x.size() < background_dimension) throw std::invalid_argument("The length of the vector argument of Geometry::set_coordinates must not be less than the background dimension!");

    unsigned int i;
    std::vector<kind> xt;
    unsigned int ulimit = (uniform) ? background_dimension : x.size();
    
    for(i=0; i<ulimit; ++i) {
      xt.push_back(invert(x[i]));
    }
    coordinates[v] = xt;
  }

  template<class kind>
  inline double Geometry<kind>::get_squared_distance(int v,const std::vector<double>& x) const 
  {
    if (relational) throw std::runtime_error("Illegal method call (Geometry::get_squared_distance) for relational model!");
    if (x.size() < background_dimension) throw std::invalid_argument("The length of the vector argument of Geometry::get_squared_distance must not be less than the background dimension!");

    unsigned int i;    
    std::vector<double> base;

    for(i=0; i<coordinates[v].size(); ++i) {
      base.push_back(convert(coordinates[v][i]));
    }

    double delta = (base[0] - x[0])*(base[0] - x[0]);
    if (!euclidean) delta = -delta;

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

  template<class kind>
  inline double Geometry<kind>::get_squared_distance(int v1,int v2,bool recompute) const
  {
    if (v1 == v2) throw std::invalid_argument("The vertex arguments in Geometry::get_squared_distance must be distinct!");
    if (relational && recompute) throw std::invalid_argument("The squared distances are fundamental in a relational geometry!");

    double l = 0.0;

    if (relational || (high_memory && !recompute)) {
      l = convert(distances[compute_index(v1,v2)]);
      if (std::isnan(l)) throw std::runtime_error("NaN detected in Geometry::get_squared_distance for vertices " + std::to_string(v1) + " and " + std::to_string(v2)); 
      return l;
    }

    unsigned int n,m;
    if (uniform) {
      m = background_dimension;
    }
    else {
      m = coordinates[v1].size();
      n = coordinates[v2].size();
      m = (m <= n) ? m : n;
    }

    kind q = (coordinates[v1][0] - coordinates[v2][0])*(coordinates[v1][0] - coordinates[v2][0]);
    if (!euclidean) q = -q;
    for(n=1; n<m; ++n) {
      q += (coordinates[v1][n] - coordinates[v2][n])*(coordinates[v1][n] - coordinates[v2][n]);
    }
    l = convert(q,2);
    return l;
  }
}
#endif


