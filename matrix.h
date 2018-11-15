#include "global.h"

#ifndef _matrixh
#define _matrixh

namespace SYNARMOSMA {
  template<class kind>
  class Matrix;

  template<class kind>
  std::ostream& operator <<(std::ostream&,const Matrix<kind>&);

  template<class kind>
  std::vector<kind>& operator *(const Matrix<kind>&,const std::vector<kind>&);

  template<class kind>
  Matrix<kind>& operator *(const Matrix<kind>&,const Matrix<kind>&);

  template<class kind>
  Matrix<kind>& operator +(const Matrix<kind>&,const Matrix<kind>&);

  template<class kind>
  void permute(Matrix<kind>&,unsigned int,unsigned int,char);

  template<class kind>
  void invert(Matrix<kind>&,unsigned int,char);

  template<class kind>
  void combine(Matrix<kind>&,unsigned int,unsigned int,const kind&,char);

  template<class kind>
  void pivot_row(Matrix<kind>&,unsigned int,unsigned int);

  template<class kind>
  void pivot_column(Matrix<kind>&,unsigned int,unsigned int);

  template<class kind>
  void move_minimum(Matrix<kind>&,unsigned int);

  template<class kind>
  void prepare_matrix(Matrix<kind>&,unsigned int);

  template<class kind>
  unsigned int normalize(Matrix<kind>&);

  template<class kind>
  /// A class representing a general rectangular matrix over an arbitrary ring or field of elements, using templates.
  class Matrix {
   protected:
    /// This property is heart of the Matrix class as it contains all of the 
    /// elements stored in compressed form as an array of vectors. There is a 
    /// vector for each row, with the vector's elements consisting of a pair. 
    /// The first part of the pair is the value of the matrix element, the 
    /// second is the column index.
    std::vector<std::pair<kind,unsigned int> >* elements;
    /// The number of rows of this  
    /// matrix.
    unsigned int nrow = 0;
    /// The number of columns of this 
    /// matrix
    unsigned int ncolumn = 0;
    /// A property which determines whether or not the matrix 
    /// has been transformed into row-reduced form by the normalize() 
    /// function.
    bool normalized = false;

    void display(std::ostream&) const;
    bool divisible(unsigned int,unsigned int*,unsigned int*,kind*) const;
    /// This method computes the transpose of the matrix provided as the argument, setting the output to the current matrix.
    void transpose(const Matrix<kind>&);
    int write_elements(std::ofstream&) const;
    int read_elements(std::ifstream&);
   public:
    static const kind zero;
    static const kind neg1;
    static const kind unity;

    /// The default constructor which does nothing in the case of this class.
    Matrix();
    /// The constructor for a square matrix whose dimension is the first argument and whose second optional argument is to construct the identity matrix when true.
    Matrix(unsigned int,bool = false);
    /// The constructor for a general rectangular matrix, which is initialized to the zero matrix.
    Matrix(unsigned int,unsigned int);
    Matrix(const Matrix<kind>&);
    /// The destructor for this class - if Matrix::nrow is greater than zero, it frees the memory in the Matrix::elements property.
    ~Matrix();
    /// This method sets Matrix::normalized to false and, if its argument is true, checks if Matrix::nrow is greater than zero and if so frees the memory in Matrix::elements, then sets Matrix::nrow and Matrix::ncolumn to zero; if the argument is false the method just clears the values of each row vector.
    void clear(bool = true);
    /// This method checks the internal consistency of the matrix: Matrix::nrow and Matrix::ncolumn must both be greater than zero, and every column value in the Matrix::elements property must lie between 0 and Matrix::ncolumn-1 and must only occur once. If any of these conditions are not satisfied, the method returns false.
    bool consistent() const;
    /// This method writes the instance properties to a binary disk file and returns the number of bytes written to the file.
    int serialize(std::ofstream&) const;
    /// This method calls the clear() method on the instance and then reads the properties from a binary disk file and returns the number of bytes read.
    int deserialize(std::ifstream&);
    /// This method sets the value of Matrix::nrow and Matrix::ncolumn to the two arguments and then allocates the memory for the Matrix::elements property.    
    void initialize(unsigned int,unsigned int);
    /// This method sets the value of the matrix element specified by the first two arguments to the third argument and, if the fourth argument is true, increments the value rather than overwriting it.
    inline void set(unsigned int,unsigned int,kind,bool = false);
    /// This method gets the value of the element specified by the two arguments, the row and column index respectively. 
    inline kind get(unsigned int,unsigned int) const;
    /// This method increments the value of the element, specified by the two first arguments, by the third argument; if the element specified does not exist, this method has the same effect as the set() method.
    inline void increment(unsigned int,unsigned int,kind);
    /// This method determines if the row whose index is given by the argument contains any elements, returning true if the row is empty and false otherwise.
    inline bool empty_row(unsigned int) const;
    /// This method returns the number of rows in this matrix.
    inline unsigned int get_nrow() const {return nrow;};
    /// This method computes the matrix's density, i.e. the number of non-zero elements divided by the total number of elements, and returns this value.
    inline double density() const {return double(number_nonzero())/double(nrow*ncolumn);};
    /// This method computes the number of non-zero elements in this matrix and returns this value.
    inline unsigned int number_nonzero() const; 
    /// This method calculates the percentage of the matrix's rows which are not diagonally dominant, i.e. for which \f$|a_{ii}| < \sum_{j=1, j\ne i}^N |a_{ij}|\f$. 
    double dispersion() const;
    void convert(kind*,char) const;
    /// This method obtains the vector of diagonal elements of the matrix, i.e. those elements whose row index is the same as their column index; the output vector will have a length of Matrix::nrow.
    void get_diagonal(std::vector<kind>&) const;
    bool diagonally_dominant() const;
    bool diagonally_dominant(unsigned int) const;
    bool diagonally_dominant(unsigned int,unsigned int) const;
    unsigned int optimize_dominance(std::vector<unsigned int>&);
    kind determinant() const;
    void multiply(const std::vector<kind>&,std::vector<kind>&) const;
    void increment(const Matrix<kind>&);
    void homotopy_scaling(kind,kind,Matrix<kind>*) const;
    void gauss_seidel_solver(std::vector<kind>&,const std::vector<kind>&,double,int);
    kind get_first_nonzero(unsigned int) const;
    void get_row(std::vector<kind>&,std::vector<unsigned int>&,unsigned int) const;
    Matrix<kind>& operator =(const Matrix<kind>&);
    friend void permute<>(Matrix<kind>&,unsigned int,unsigned int,char);
    friend void invert<>(Matrix<kind>&,unsigned int,char);
    friend void combine<>(Matrix<kind>&,unsigned int,unsigned int,const kind&,char);
    friend void pivot_row<>(Matrix<kind>&,unsigned int,unsigned int);
    friend void pivot_column<>(Matrix<kind>&,unsigned int,unsigned int);
    friend void move_minimum<>(Matrix<kind>&,unsigned int);
    friend void prepare_matrix<>(Matrix<kind>&,unsigned int);
    friend unsigned int normalize<>(Matrix<kind>&);
    friend std::ostream& operator << <>(std::ostream&,const Matrix<kind>&);
  };

  template<class kind>
  std::ostream& operator <<(std::ostream& os,const Matrix<kind>& source)
  {
    unsigned int i,j,k;
    kind w;
    bool found;

    for(i=0; i<source.nrow; ++i) {
      for(j=0; j<source.ncolumn; ++j) {
        w = 0;
        found = false;
        for(k=0; k<source.elements[i].size(); ++k) {
          if (source.elements[i][k].second == j) {
            w = source.elements[i][k].first;
            found = true;
            break;
          }
        }
        if (!found) {
          os << "%" << w << " ";
        }
        else {
          os << w << " ";
        }
      }
      if (i < (source.nrow - 1)) os << std::endl;
    }
    return os;
  }

  template<class kind>
  kind Matrix<kind>::get(unsigned int n,unsigned int m) const
  {
    if (n >= nrow) throw std::invalid_argument("The row number argument is illegal for this matrix!");
    if (m >= ncolumn) throw std::invalid_argument("The column number argument is illegal for this matrix!");
    unsigned int i,mu = m;
    kind output = zero;

    if (elements[n].empty()) return output;
    for(i=0; i<elements[n].size(); ++i) {
      if (elements[n][i].second == mu) {
        output = elements[n][i].first;
        break;
      } 
    }
    return output;
  }

  template<class kind>
  void Matrix<kind>::set(unsigned int n,unsigned int m,kind v,bool increment)
  {
    if (n >= nrow) throw std::invalid_argument("The row number argument is illegal for this matrix!");
    if (m >= ncolumn) throw std::invalid_argument("The column number argument is illegal for this matrix!");
    unsigned int i,q;
    bool found = false;
    for(i=0; i<elements[n].size(); ++i) {
      if (elements[n][i].second == m) {
        found = true;
        q = i;
        break;
      }
    }
    if (found) {
      if (increment) {
        elements[n][q].first += v;
      }
      else {
        elements[n][q].first = v;
      }
    }
    else {
      elements[n].push_back(std::pair<kind,unsigned int>(v,m));
    }
  }

  template<class kind>
  void Matrix<kind>::increment(unsigned int n,unsigned int m,kind v)
  {
    if (n >= nrow) throw std::invalid_argument("The row number argument is illegal for this matrix!");
    if (m >= ncolumn) throw std::invalid_argument("The column number argument is illegal for this matrix!");
    unsigned int i,q,mu = m;
    bool found = false;

    for(i=0; i<elements[n].size(); ++i) {
      if (elements[n][i].second == mu) {
        found = true;
        q = i;
        break;
      }
    }
    if (found) {
      elements[n][q].first += v;
    }
    else {
      elements[n].push_back(std::pair<kind,unsigned int>(v,m));
    }
  }

  template<class kind>
  unsigned int Matrix<kind>::number_nonzero() const
  {
    unsigned int i,nzero = 0;

    for(i=0; i<nrow; ++i) {
      nzero += elements[i].size();
    }
    return nzero;
  }

  template<class kind>
  bool Matrix<kind>::empty_row(unsigned int r) const
  {
    if (r >= nrow) throw std::invalid_argument("The row number argument is illegal for this matrix!");
    return elements[r].empty();
  }
}
#endif
