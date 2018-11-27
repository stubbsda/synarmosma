#include "random.h"

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
    /// has been transformed into reduced row echelon form by the normalize() 
    /// method.
    bool normalized = false;

    /// This method writes out the complete matrix to the output stream, one line per row, with a '%' before the entries explicitly stored in the Matrix::elements vector for that row.
    void display(std::ostream&) const;
    /// This method iterates through all the putatively non-zero matrix elements and if it finds one which is indistinguishable from zero, it removes this element; the method returns the number of such pseudo-elements found in the matrix.
    unsigned int eliminate_zeros();
    /// This method checks whether the diagonal element (if zero the method returns false) of the row whose index is specified by the first argument is divisible by an element of the matrix with a row and column index greater than it, returning true in that case along with the indices of the element and the quotient. For a base type like double or std::complex<double> the divisibility test is of course meaningless.
    bool divisible(unsigned int,unsigned int*,unsigned int*,kind*) const;
    /// This method computes the transpose of the matrix provided as the argument, setting the output to the current matrix.
    void transpose(const Matrix<kind>&);
    /// This method writes the content of the matrix itself in binary format to an output stream and is used by the serialize() method; it returns the number of bytes written.
    int write_elements(std::ofstream&) const;
    /// This method reads the content of the matrix itself in binary format from an input stream is used by the deserialize() method; it returns the number of bytes read.
    int read_elements(std::ifstream&);
   public:
    /// The value 0, stored in the correct data type for this instantiation of 
    /// the template class.
    static const kind zero;
    /// The value -1, stored in the correct data type for this instantiation of 
    /// the template class.
    static const kind neg1;
    /// The value 1, stored in the correct data type for this instantiation of 
    /// the template class.
    static const kind unity;

    /// The default constructor which does nothing in the case of this class.
    Matrix();
    /// The constructor for a square matrix whose dimension is the first argument and whose second optional argument is to construct the identity matrix when true.
    Matrix(unsigned int,bool = false);
    /// The constructor for a general rectangular matrix, which is initialized to the zero matrix.
    Matrix(unsigned int,unsigned int);
    /// The standard copy constructor - it calls the clear() method and then copies over all properties from the soure matrix to this one.
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
    /// This method returns the number of rows in this matrix.
    inline unsigned int get_nrow() const {return nrow;};
    /// This method computes the matrix's density, i.e. the number of non-zero elements divided by the total number of elements, and returns this value.
    inline double density() const {return double(number_nonzero())/double(nrow*ncolumn);};
    /// This method computes the number of non-zero elements in this matrix and returns this value.
    unsigned int number_nonzero() const; 
    /// This method determines if the row whose index is given by the argument contains any elements, returning true if the row is empty and false otherwise.
    bool empty_row(unsigned int) const;
    /// This method gets the non-zero element of the row with the lowest column index and returns it; if the row has no non-zero elements, it returns zero. 
    kind get_first_nonzero(unsigned int) const;
    /// This method calculates the percentage of the matrix's rows which are not diagonally dominant, i.e. for which \f$|a_{ii}| < \sum_{j=1, j\ne i}^N |a_{ij}|\f$. 
    double dispersion() const;
    /// This method writes the matrix into a one-dimensional array (the first argument) whose length is the product of Matrix::nrow and Matrix::ncolumn, using either a row-oriented or column-oriented convention, specified by the second argument ('r' or 'c').
    void convert(kind*,char) const;
    /// This method extracts the diagonal element of the row whose index is the argument.
    inline kind get_diagonal(unsigned int) const;
    /// This method obtains the vector of diagonal elements of the matrix, i.e. those elements whose row index is the same as their column index; the output vector will have a length of Matrix::nrow.
    void get_diagonal(std::vector<kind>&) const;
    /// This method checks if the matrix is diagonally dominant, i.e. if for every row \f$i\f$ the inequality \f$|a_{ii}| \ge \sum_{j=1, j\ne i}^N |a_{ij}|\f$ is satisfied, and returns true if this is the case.
    bool diagonally_dominant() const;
    /// This method checks if the row whose index is this method's unique argument is diagonally dominant, i.e. \f$|a_{ii}| \ge \sum_{j=1, j\ne i}^N |a_{ij}|\f$ for the row \f$i\f$, returning true if this is the case.
    bool diagonally_dominant(unsigned int) const;
    /// This method checks if the element specified by the two arguments (row and column index) dominates its row, i.e. \f$|a_{ik}| \ge \sum_{j=1, j\ne k}^N |a_{ij}|\f$, where \f$i\f$ and \f$k\f$ are the two arguments of this method. 
    bool diagonally_dominant(unsigned int,unsigned int) const;
    /// This method permutes the rows of the matrix in an effort to maximize the number of rows which are diagonally dominant; the argument is filled with the row permutations and the output is the number of rows which are diagonally dominant.
    unsigned int optimize_dominance(std::vector<unsigned int>&);
    /// This method computes and returns the determinant of the matrix, first checking that it is square and then calculating it by Laplace's formula, i.e. the sum of cofactors. 
    kind determinant() const;
    /// This method multiplies the matrix by the first argument of this method and writes the output into the second argument, after checking that the vector conforms to the matrix dimensions.
    void multiply(const std::vector<kind>&,std::vector<kind>&) const;
    /// This method adds the matrix that is the method's argument to the instance, after first checking that the dimensions of the two matrices are identical.
    void increment(const Matrix<kind>&);
    /// This method multiplies the instance matrix by the first argument and then adds to it the identity matrix of appropriate dimension multiplied by the second argument, writing the result in the final argument; in summary, \f$B = \lambda_1 A + \lambda_2 I\f$ where the three arguments are \f$\lambda_1\f$, \f$\lambda_2\f$ and \f$B\f$, with \f$A\f$ the current instance.
    void homotopy_scaling(kind,kind,Matrix<kind>*) const;
    /// This method uses the Gauss-Seidel iterative algorithm to attempt to solve the linear system \f$A x = b\f$, where \f$A\f$ is the current instance, \f$x\f$ and \f$b\f$ the first and second arguments. The third argument is the convergence threshold and the fourth the maximum number of iterations to perform. 
    void gauss_seidel_solver(std::vector<kind>&,const std::vector<kind>&,double,int);
    /// This method takes as its input the final argument, a row index, and then puts the content of this row in the two initial arguments, the first containing the values and the second the column indices.
    void get_row(std::vector<kind>&,std::vector<unsigned int>&,unsigned int) const;
    /// The overloaded assignment operator for this class, which first calls clear() and then behaves exactly like the copy constructor for this class.
    Matrix<kind>& operator =(const Matrix<kind>&);
    /// This global function permutes two distinct columns (when the final argument is 'c') or two distinct rows (when this argument is 'r'), whose indices are the second and third arguments; the matrix affected is the first argument.
    friend void permute<>(Matrix<kind>&,unsigned int,unsigned int,char);
    /// This global function multiplies the elements of a row (when the final argument is 'r') or column (when it is 'c') by -1, with the index specified by the second argument and the matrix itself the first argument.
    friend void invert<>(Matrix<kind>&,unsigned int,char);
    /// This global function performs a scaled addition of two rows (when the final argument is 'r') or two columns (when it is 'c') of the matrix that is the first argument. As an example when the final argument is 'r', this method performs \f$A_{nk} = A_{nk} + q A_{mk}\f$ for all \f$k\f$, where \f$n\f$ and \f$m\f$ are the second and third arguments, with the scalar \f$q\f$ the fourth argument.
    friend void combine<>(Matrix<kind>&,unsigned int,unsigned int,const kind&,char);
    /// This global function accepts a row and column index (the final two arguments) and uses the combine() function to add two rows to eliminate an off-diagonal element.
    friend void pivot_row<>(Matrix<kind>&,unsigned int,unsigned int);
    /// This global function accepts a row and column index (the final two arguments) and uses the combine() function to add two columns to eliminate an off-diagonal element.
    friend void pivot_column<>(Matrix<kind>&,unsigned int,unsigned int);
    /// This global function uses the permute() function to move the row whose index is the second argument so that the smallest non-zero element in the sub-matrix with row and column index greater than or equal to the second argument.
    friend void move_minimum<>(Matrix<kind>&,unsigned int);
    /// This global function works on a particular row (whose index is the second argument) of a matrix, calling the combine(), move_minimum(), pivot_row() and pivot_column() functions to carry out elementary row operations.
    friend void prepare_matrix<>(Matrix<kind>&,unsigned int);
    /// This global function transforms its unique argument into reduced row echelon form by elementary row operations, calling the prepare_matrix() and invert() functions for this purpose. The return value is the number of diagonal entries equal to 1.
    friend unsigned int normalize<>(Matrix<kind>&);
    /// This overloading of the ostream operator writes the matrix to the screen according to the same conventions as the display() method.
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
    kind output = Matrix<kind>::zero; 

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
  kind Matrix<kind>::get_diagonal(unsigned int n) const
  {
    if (n >= nrow) throw std::invalid_argument("The row number argument is illegal for this matrix!");
    unsigned int i;
    kind output = Matrix<kind>::zero;

    for(i=0; i<elements[n].size(); ++i) {
      if (elements[n][i].second == n) {
        output = elements[n][i].first;
        break;
      }
    }
    return output;
  }
}
#endif
