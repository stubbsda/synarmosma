#include "global.h"

#ifndef _matrixh
#define _matrixh

namespace SYNARMOSMA {
  template<class kind>
  class Matrix;

  template<class kind>
  std::ostream& operator <<(std::ostream&,const Matrix<kind>&);

  template<class kind>
  Matrix<kind> operator -(const Matrix<kind>&);

  template<class kind>
  Matrix<kind> operator +(const Matrix<kind>&,const Matrix<kind>&);

  template<class kind>
  Matrix<kind> operator -(const Matrix<kind>&,const Matrix<kind>&);

  template<class kind>
  Matrix<kind> operator *(kind,const Matrix<kind>&);

  template<class kind>
  std::vector<kind> operator *(const Matrix<kind>&,const std::vector<kind>&);

  template<class kind>
  Matrix<kind> operator *(const Matrix<kind>&,const Matrix<kind>&);

  template<class kind>
  Matrix<kind> operator ^(const Matrix<kind>&,int);

  /// A template class representing a general rectangular matrix over a floating point base type.
  template<class kind>
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

    /// This method writes out the complete matrix to the output stream, one line per row, with a '%' before the entries explicitly stored in the Matrix::elements vector for that row.
    void display(std::ostream&) const;
    /// This method iterates through all the putatively non-zero matrix elements and if it finds one which is indistinguishable from zero, it removes this element; the method returns the number of such pseudo-elements found in the matrix.
    unsigned int eliminate_zeros();
    /// This method computes the transpose of the matrix provided as the argument, setting the output to the current matrix.
    void transpose(const Matrix<kind>&);
    /// This method writes the content of the matrix itself in binary format to an output stream and is used by the serialize() method; it returns the number of bytes written.
    int write_elements(std::ofstream&) const;
    /// This method reads the content of the matrix itself in binary format from an input stream is used by the deserialize() method; it returns the number of bytes read.
    int read_elements(std::ifstream&);
   public:
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
    void set(unsigned int,unsigned int,kind,bool = false);
    /// This method gets the value of the element specified by the two arguments, the row and column index respectively. 
    kind get(unsigned int,unsigned int) const;
    /// This method increments the value of the element, specified by the two first arguments, by the third argument; if the element specified does not exist, this method has the same effect as the set() method.
    void increment(unsigned int,unsigned int,kind);
    /// This method returns the number of rows in this matrix.
    unsigned int get_nrow() const;
    /// This method returns the number of rows in this matrix.
    unsigned int get_ncolumn() const;
    /// This method computes the matrix's density, i.e. the number of non-zero elements divided by the total number of elements, and returns this value.
    double density() const;
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
    kind get_diagonal(unsigned int) const;
    /// This method obtains the vector of diagonal elements of the matrix, i.e. those elements whose row index is the same as their column index; the output vector will have a length of Matrix::nrow.
    void get_diagonal(std::vector<kind>&) const;
    /// This method checks if the matrix is diagonally dominant, i.e. if for every row \f$i\f$ the inequality \f$|a_{ii}| \ge \sum_{j=1, j\ne i}^N |a_{ij}|\f$ is satisfied, and returns true if this is the case. The method's unique argument controls whether or not any output concerning the matrix's diagonal dominance is written to the console.
    bool diagonally_dominant(bool = false) const;
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
    /// This method uses the Gauss-Seidel iterative algorithm to attempt to solve the linear system \f$A x = b\f$, where \f$A\f$ is the current instance, \f$x\f$ and \f$b\f$ the first and second arguments. The third argument is the convergence threshold and the fourth the maximum number of iterations to perform. The final argument controls whether or not any output concerning the progress of the iterative solver is written to the console.
    void gauss_seidel_solver(std::vector<kind>&,const std::vector<kind>&,double,int,bool = false);
    /// This method takes as its input the final argument, a column index, and then writes the complete column vector (zeros included) in the first argument.
    void get_column(std::vector<kind>&,unsigned int) const;
    /// This method takes as its input the final argument, a row index, and then writes the complete row vector (zeros included) in the first argument.
    void get_row(std::vector<kind>&,unsigned int) const;
    /// This method takes as its input the final argument, a row index, and then puts the content of this row in the two initial arguments, the first containing the values and the second the column indices.
    void get_row(std::vector<kind>&,std::vector<unsigned int>&,unsigned int) const;
    /// This method takes as its input the final argument, a row index, and then changes the current values for this row to the non-zero elements of the first argument, which must be of length Matrix::nrow.
    void set_row(const std::vector<kind>&,unsigned int);
    /// This method computes and returns the value of the matrix raised to the power that is the method's argument, i.e. \f$B = A^n\f$.
    Matrix<kind> pow(unsigned int) const;
    /// The overloaded assignment operator for this class, which first calls clear() and then behaves exactly like the copy constructor for this class.
    Matrix<kind>& operator =(const Matrix<kind>&);
    /// This overloading of the ostream operator writes the matrix to the screen according to the same conventions as the display() method.
    friend std::ostream& operator << <>(std::ostream&,const Matrix<kind>&);
    /// This overloading of the multiplication operator carries out the usual multiplication of a matrix by a scalar, adapted to the compressed storage format used by this class.
    friend Matrix<kind> operator *<>(kind,const Matrix<kind>&);
    /// This overloading of the multiplication operator carries out the usual matrix multiplication, adapted to the compressed storage format used by this class.
    friend Matrix<kind> operator *<>(const Matrix<kind>&,const Matrix<kind>&);
    /// This overloading of the multiplication operator carries out the usual matrix-vector multiplication, adapted to the compressed storage format used by this class.
    friend std::vector<kind> operator *<>(const Matrix<kind>&,const std::vector<kind>&);
    /// This overloading of the addition operator carries out the usual element-wise matrix addition, adapted to the compressed storage format used by this class.
    friend Matrix<kind> operator +<>(const Matrix<kind>&,const Matrix<kind>&);
  };

  template<class kind>
  Matrix<kind> operator *(kind alpha,const Matrix<kind>& A)
  {
    unsigned int i,j;
    Matrix<kind> output = A;

    for(i=0; i<output.nrow; ++i) {
      for(j=0; j<output.elements[i].size(); ++j) {
        output.elements[i][j].first *= alpha;
      }      
    }
    output.eliminate_zeros();
    return output;
  }

  template<class kind>
  std::vector<kind> operator *(const Matrix<kind>& A,const std::vector<kind>& b)
  {
    if (A.ncolumn != b.size()) throw std::invalid_argument("The matrix and vector dimensions don't conform for multiplication!");
    unsigned int i,j;
    kind sum;
    std::vector<kind> output;

    for(i=0; i<A.nrow; ++i) {
      sum = kind(0.0);
      for(j=0; j<A.elements[i].size(); ++j) {
        sum += A.elements[i][j].first*b[A.elements[i][j].second];
      }
      output.push_back(sum);
    }

    return output;
  }

  template<class kind>
  Matrix<kind> operator -(const Matrix<kind>& p)
  {
    Matrix<kind> output = kind(-1)*p;
    
    return output;
  }

  template<class kind>
  Matrix<kind> operator -(const Matrix<kind>& p1,const Matrix<kind>& p2)
  {
    Matrix<kind> output = p1 + (-p2);
    
    return output;
  }


  template<class kind>
  Matrix<kind> operator +(const Matrix<kind>& A,const Matrix<kind>& B)
  {
    if (A.nrow != B.nrow || A.ncolumn != B.ncolumn) throw std::invalid_argument("The dimensions of the two matrices don't conform for addition!");
    unsigned int i,j,k,in1,cvalue;
    bool found;
    Matrix<kind> output(A.nrow,A.ncolumn);

    output = A;

    // Now the matrix multiplication itself...
    for(i=0; i<B.nrow; ++i) {
      for(j=0; j<B.elements[i].size(); ++j) {
        in1 = B.elements[i][j].second;
        found = false;
        for(k=0; k<output.elements[i].size(); ++k) {
          if (output.elements[i][k].second == in1) {
            cvalue = k;
            found = true;
            break;
          }
        }
        if (found) {
          output.elements[i][cvalue].first += B.elements[i][j].first;        
        }
        else {
          output.elements[i].push_back(B.elements[i][j]);
        }
      }
    }
    output.eliminate_zeros();
    return output;
  }

  template<class kind>
  Matrix<kind> operator *(const Matrix<kind>& A,const Matrix<kind>& B)
  {
    if (A.ncolumn != B.nrow) throw std::invalid_argument("The dimensions of the two matrices don't conform for multiplication!");
    unsigned int i,j,k,l,in1;
    kind sum;
    Matrix<kind> output(A.nrow,B.ncolumn);

    // Convert the matrix B to a column-oriented format...
    std::vector<kind>* BC = new std::vector<kind>[B.ncolumn];
    std::vector<unsigned int>* BC_row = new std::vector<unsigned int>[B.ncolumn];
    for(l=0; l<B.ncolumn; ++l) {
      for(i=0; i<B.nrow; ++i) {
        for(j=0; j<B.elements[i].size(); ++j) {
          if (B.elements[i][j].second == l) {
            BC[l].push_back(B.elements[i][j].first);
            BC_row[l].push_back(i);
          }
        }
      }
    }
    // Now the matrix multiplication itself...
    for(i=0; i<A.nrow; ++i) {
      for(j=0; j<B.ncolumn; ++j) {
        sum = kind(0.0);
        for(k=0; k<A.elements[i].size(); ++k) {
          in1 = A.elements[i][k].second;
          for(l=0; l<BC_row[j].size(); ++l) {
            if (BC_row[j][l] == in1) {
              sum += A.elements[i][k].first*BC[j][l];
            }
          }
        }
        output.elements[i].push_back(std::pair<kind,unsigned int>(sum,j));
      }
    }
    output.eliminate_zeros();
    delete[] BC;
    delete[] BC_row;
    return output;
  }

  template<class kind>
  Matrix<kind> operator ^(const Matrix<kind>& p,int n)
  {
    if (n < 1) throw std::invalid_argument("The exponent must be positive!");

    Matrix<kind> output = p;

    for(int i=1; i<n; ++i) {
      output = output*p;
    }

    return output;
  }

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
  inline unsigned int Matrix<kind>::get_nrow() const
  {
    return nrow;
  }

  template<class kind>
  inline unsigned int Matrix<kind>::get_ncolumn() const
  {
    return ncolumn;
  }

  template<class kind>
  inline double Matrix<kind>::density() const
  {
    return double(number_nonzero())/double(nrow*ncolumn);
  }

  template<class kind>
  inline kind Matrix<kind>::get(unsigned int n,unsigned int m) const
  {
    if (n >= nrow) throw std::invalid_argument("The row number argument is illegal for this matrix!");
    if (m >= ncolumn) throw std::invalid_argument("The column number argument is illegal for this matrix!");
    unsigned int i,mu = m;
    kind output = kind(0.0); 

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
  inline void Matrix<kind>::set(unsigned int n,unsigned int m,kind v,bool increment)
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
  inline void Matrix<kind>::increment(unsigned int n,unsigned int m,kind v)
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
  inline kind Matrix<kind>::get_diagonal(unsigned int n) const
  {
    if (n >= nrow) throw std::invalid_argument("The row number argument is illegal for this matrix!");
    unsigned int i;
    kind output = kind(0.0);

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
