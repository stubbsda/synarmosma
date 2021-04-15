#include "global.h"

#ifndef _integermatrixh
#define _integermatrixh

namespace SYNARMOSMA {
  template<class kind>
  class Integer_Matrix;

  template<class kind>
  std::ostream& operator <<(std::ostream&,const Integer_Matrix<kind>&);

  template<class kind>
  std::vector<kind> operator *(const Integer_Matrix<kind>&,const std::vector<kind>&);

  template<class kind>
  Integer_Matrix<kind> operator *(const Integer_Matrix<kind>&,const Integer_Matrix<kind>&);

  template<class kind>
  Integer_Matrix<kind> operator +(const Integer_Matrix<kind>&,const Integer_Matrix<kind>&);

  template<class kind>
  void permute(Integer_Matrix<kind>&,unsigned int,unsigned int,char);

  template<class kind>
  void invert(Integer_Matrix<kind>&,unsigned int,char);

  template<class kind>
  void combine(Integer_Matrix<kind>&,unsigned int,unsigned int,const kind&,char);

  template<class kind>
  void pivot_row(Integer_Matrix<kind>&,unsigned int,unsigned int);

  template<class kind>
  void pivot_column(Integer_Matrix<kind>&,unsigned int,unsigned int);

  template<class kind>
  void move_minimum(Integer_Matrix<kind>&,unsigned int);

  template<class kind>
  void prepare_matrix(Integer_Matrix<kind>&,unsigned int);

  template<class kind>
  unsigned int normalize(Integer_Matrix<kind>&);

  /// A template class representing a general rectangular matrix over an integer base type.
  template<class kind>
  class Integer_Matrix {
   protected:
    /// This property is heart of the Integer_Matrix class as it contains all of the
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
    /// A property which determines whether or not the matrix has been transformed
    /// into reduced row echelon form by the normalize() method.
    bool normalized = false;

    /// This method writes out the complete matrix to the output stream, one line per row, with a '%' before the entries explicitly stored in the Integer_Matrix::elements vector for that row.
    void display(std::ostream&) const;
    /// This method iterates through all the putatively non-zero matrix elements and if it finds one which is indistinguishable from zero, it removes this element; the method returns the number of such pseudo-elements found in the matrix.
    unsigned int eliminate_zeros();
    /// This method checks whether the diagonal element (if zero the method returns false) of the row whose index is specified by the first argument is divisible by an element of the matrix with a row and column index greater than it, returning true in that case along with the indices of the element and the quotient. For a base type like double or std::complex<double> the divisibility test is of course meaningless.
    bool divisible(unsigned int,unsigned int*,unsigned int*,kind*) const;
    /// This method computes the transpose of the matrix provided as the argument, setting the output to the current matrix.
    void transpose(const Integer_Matrix<kind>&);
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
    Integer_Matrix();
    /// The constructor for a square matrix whose dimension is the first argument and whose second optional argument is to construct the identity matrix when true.
    Integer_Matrix(unsigned int,bool = false);
    /// The constructor for a general rectangular matrix, which is initialized to the zero matrix.
    Integer_Matrix(unsigned int,unsigned int);
    /// The standard copy constructor - it calls the clear() method and then copies over all properties from the source matrix to this one.
    Integer_Matrix(const Integer_Matrix<kind>&);
    /// The destructor for this class - if Integer_Matrix::nrow is greater than zero, it frees the memory in the Integer_Matrix::elements property.
    ~Integer_Matrix();
    /// This method sets Integer_Matrix::normalized to false and, if its argument is true, checks if Integer_Matrix::nrow is greater than zero and if so frees the memory in Integer_Matrix::elements, then sets Integer_Matrix::nrow and Integer_Matrix::ncolumn to zero; if the argument is false the method just clears the values of each row vector.
    void clear(bool = true);
    /// This method checks the internal consistency of the matrix: Integer_Matrix::nrow and Integer_Matrix::ncolumn must both be greater than zero, and every column value in the Integer_Matrix::elements property must lie between 0 and Integer_Matrix::ncolumn-1 and must only occur once. If any of these conditions are not satisfied, the method returns false.
    bool consistent() const;
    /// This method writes the instance properties to a binary disk file and returns the number of bytes written to the file.
    int serialize(std::ofstream&) const;
    /// This method calls the clear() method on the instance and then reads the properties from a binary disk file and returns the number of bytes read.
    int deserialize(std::ifstream&);
    /// This method sets the value of Integer_Matrix::nrow and Integer_Matrix::ncolumn to the two arguments and then allocates the memory for the Integer_Matrix::elements property.    
    void initialize(unsigned int,unsigned int);
    /// This method sets the value of the matrix element specified by the first two arguments to the third argument and, if the fourth argument is true, increments the value rather than overwriting it.
    void set(unsigned int,unsigned int,kind,bool = false);
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
    /// This method gets the value of the element specified by the two arguments, the row and column index respectively. 
    kind get(unsigned int,unsigned int) const;
    /// This method extracts the diagonal element of the row whose index is the argument.
    kind get_diagonal(unsigned int) const;
    /// This method obtains the vector of diagonal elements of the matrix, i.e. those elements whose row index is the same as their column index; the output vector will have a length of Integer_Matrix::nrow.
    void get_diagonal(std::vector<kind>&) const;
    /// This method computes whether or not this instance of the class is a symmetric matrix (i.e. it is square and additionally equal to its transpose), returning true if this is so and false otherwise.
    bool symmetric() const;
    /// This method computes and returns the determinant of the matrix, first checking that it is square and then calculating it by Laplace's formula, i.e. the sum of cofactors. 
    kind determinant() const;
    /// This method multiplies the matrix by the first argument of this method and writes the output into the second argument, after checking that the vector conforms to the matrix dimensions.
    void multiply(const std::vector<kind>&,std::vector<kind>&) const;
    /// This method adds the matrix that is the method's argument to the instance, after first checking that the dimensions of the two matrices are identical.
    void increment(const Integer_Matrix<kind>&);
    /// This method takes as its input the final argument, a row index, and then puts the content of this row in the two initial arguments, the first containing the values and the second the column indices.
    void get_row(std::vector<kind>&,std::vector<unsigned int>&,unsigned int) const;
    /// The overloaded assignment operator for this class, which first calls clear() and then behaves exactly like the copy constructor for this class.
    Integer_Matrix<kind>& operator =(const Integer_Matrix<kind>&);
    /// This global function permutes two distinct columns (when the final argument is 'c') or two distinct rows (when this argument is 'r'), whose indices are the second and third arguments; the matrix affected is the first argument.
    friend void permute<>(Integer_Matrix<kind>&,unsigned int,unsigned int,char);
    /// This global function multiplies the elements of a row (when the final argument is 'r') or column (when it is 'c') by -1, with the index specified by the second argument and the matrix itself the first argument.
    friend void invert<>(Integer_Matrix<kind>&,unsigned int,char);
    /// This global function performs a scaled addition of two rows (when the final argument is 'r') or two columns (when it is 'c') of the matrix that is the first argument. As an example when the final argument is 'r', this method performs \f$A_{nk} = A_{nk} + q A_{mk}\f$ for all \f$k\f$, where \f$n\f$ and \f$m\f$ are the second and third arguments, with the scalar \f$q\f$ the fourth argument.
    friend void combine<>(Integer_Matrix<kind>&,unsigned int,unsigned int,const kind&,char);
    /// This global function accepts a row and column index (the final two arguments) and uses the combine() function to add two rows to eliminate an off-diagonal element.
    friend void pivot_row<>(Integer_Matrix<kind>&,unsigned int,unsigned int);
    /// This global function accepts a row and column index (the final two arguments) and uses the combine() function to add two columns to eliminate an off-diagonal element.
    friend void pivot_column<>(Integer_Matrix<kind>&,unsigned int,unsigned int);
    /// This global function uses the permute() function to move the row whose index is the second argument so that the smallest non-zero element in the sub-matrix with row and column index greater than or equal to the second argument.
    friend void move_minimum<>(Integer_Matrix<kind>&,unsigned int);
    /// This global function works on a particular row (whose index is the second argument) of a matrix, calling the combine(), move_minimum(), pivot_row() and pivot_column() functions to carry out elementary row operations.
    friend void prepare_matrix<>(Integer_Matrix<kind>&,unsigned int);
    /// This global function transforms its unique argument into reduced row echelon form by elementary row operations, calling the prepare_matrix() and invert() functions for this purpose. The return value is the number of diagonal entries equal to 1.
    friend unsigned int normalize<>(Integer_Matrix<kind>&);
    /// This overloading of the ostream operator writes the matrix to the screen according to the same conventions as the display() method.
    friend std::ostream& operator << <>(std::ostream&,const Integer_Matrix<kind>&);
    /// This overloading of the multiplication operator carries out the usual matrix multiplication adapted to the compressed storage format used by this class.
    friend Integer_Matrix<kind> operator * <>(const Integer_Matrix<kind>&,const Integer_Matrix<kind>&);
    /// This overloading of the multiplication operator carries out the usual matrix-vector multiplication adapted to the compressed storage format used by this class.
    friend std::vector<kind> operator * <>(const Integer_Matrix<kind>&,const std::vector<kind>&);
    /// This overloading of the addition operator carries out the usual element-wise matrix addition adapted to the compressed storage format used by this class.
    friend Integer_Matrix<kind> operator + <>(const Integer_Matrix<kind>&,const Integer_Matrix<kind>&);
  };

  template<class kind>
  std::vector<kind> operator *(const Integer_Matrix<kind>& A,const std::vector<kind>& b)
  {
    if (A.ncolumn != b.size()) throw std::invalid_argument("The matrix and vector dimensions don't conform for multiplication!");
    unsigned int i,j;
    kind sum;
    std::vector<kind> output;

    for(i=0; i<A.nrow; ++i) {
      sum = Integer_Matrix<kind>::zero;
      for(j=0; j<A.elements[i].size(); ++j) {
        sum += A.elements[i][j].first*b[A.elements[i][j].second];
      }
      output.push_back(sum);
    }

    return output;
  }

  template<class kind>
  Integer_Matrix<kind> operator +(const Integer_Matrix<kind>& A,const Integer_Matrix<kind>& B)
  {
    if (A.nrow != B.nrow || A.ncolumn != B.ncolumn) throw std::invalid_argument("The dimensions of the two matrices don't conform for addition!");
    unsigned int i,j,k,in1,cvalue;
    bool found;
    Integer_Matrix<kind> output(A.nrow,A.ncolumn);

    output = A;

    // Now the matrix addition itself...
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
  Integer_Matrix<kind> operator *(const Integer_Matrix<kind>& A,const Integer_Matrix<kind>& B)
  {
    if (A.ncolumn != B.nrow) throw std::invalid_argument("The dimensions of the two matrices don't conform for multiplication!");
    unsigned int i,j,k,l,in1;
    kind sum;
    Integer_Matrix<kind> output(A.nrow,B.ncolumn);

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
        sum = Integer_Matrix<kind>::zero;
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
  std::ostream& operator <<(std::ostream& os,const Integer_Matrix<kind>& source)
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
  inline unsigned int Integer_Matrix<kind>::get_nrow() const
  {
    return nrow;
  }

  template<class kind>
  inline unsigned int Integer_Matrix<kind>::get_ncolumn() const
  {
    return ncolumn;
  }

  template<class kind>
  inline double Integer_Matrix<kind>::density() const
  {
    return double(number_nonzero())/double(nrow*ncolumn);
  }

  template<class kind>
  inline void Integer_Matrix<kind>::set(unsigned int n,unsigned int m,kind v,bool increment)
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
  inline void Integer_Matrix<kind>::increment(unsigned int n,unsigned int m,kind v)
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
}
#endif
