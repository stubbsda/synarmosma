#include "global.h"

#ifndef _bmatrixh
#define _bmatrixh

namespace SYNARMOSMA {
  /// A class representing a matrix whose elements are zero and one, so an arbitrary n x m matrix over the Galois field GF(2).
  class Binary_Matrix {
   private:
    /// An array of sets of unsigned integers, one for each row, which contains the index of the columns 
    /// where this row has a "1". 
    std::set<unsigned int>* elements;
    /// The number of rows of this binary 
    /// matrix.
    unsigned int nrow = 0;
    /// The number of columns of this binary 
    /// matrix.
    unsigned int ncolumn = 0;

    /// This method checks if the number of rows is greater than zero, if so it frees the memory associated with the Binary_Matrix::elements array and then sets the number of rows and columns to zero.
    void clear();
   public:
    /// The default constructor, which does nothing (i.e. no memory allocation for the array Binary_Matrix::elements).
    Binary_Matrix();
    /// The constructor for a square binary matrix, with the first argument equal to the matrix's dimension; if the second argument is false the matrix is zero, otherwise it constructs the identity matrix. 
    Binary_Matrix(unsigned int,bool = false);
    /// The constructor for a rectangular binary matrix, with the first two arguments the number of rows and columns while the optional third argument is the percentage of the (randomly chosen) matrix elements to be given a value of 1.
    Binary_Matrix(unsigned int,unsigned int,float = 0.0);
    /// Standard copy constructor that first calls the initialize() method with the dimensions of the source matrix and then copies over all the sets of column indices.  
    Binary_Matrix(const Binary_Matrix&);
    /// The overloaded assignment operator for this class, which first calls clear() and then behaves exactly like the copy constructor for this class.
    Binary_Matrix& operator =(const Binary_Matrix&);
    /// Destructor that releases the memory in the Binary_Matrix::elements property.
    ~Binary_Matrix();
    /// This method writes the instance properties to a binary disk file and returns the number of bytes written to the file.
    int serialize(std::ofstream&) const;
    /// This method calls the clear() method on the instance and then reads the properties from a binary disk file and returns the number of bytes read.
    int deserialize(std::ifstream&);
    /// This method sets the value of Binary_Matrix::nrow and Binary_Matrix::ncolumn to the two arguments and then allocates the memory for the Binary_Matrix::elements property.
    void initialize(unsigned int,unsigned int);
    /// This method computes and returns the rank of the matrix, i.e. the dimension of the vector space over GF(2) of the matrix's columns. 
    int rank() const;
    /// This method returns true if the matrix is symmetric (i.e. equal to its transpose) and false if it is non-square or non-symmetric.
    bool symmetric() const;
    /// This method returns the percentage (between 0.0 and 1.0) of matrix elements which are equal to 1.
    double density() const;
    /// This method writes the complete k-th row vector of the matrix to the second argument (assumed to be of size Binary_Matrix::ncolumn or more), where k is the first argument.
    void get_row(unsigned int,bool*) const;
    /// This method gets the indicated matrix element (row and column number), returning true if the element is 1 and false otherwise.
    inline bool get(unsigned int,unsigned int) const;
    /// This method sets the indicated matrix element (row and column number) to 1, returning true if the element had been 0, false otherwise. 
    inline bool set(unsigned int,unsigned int);
    /// This method sets the indicated matrix element (row and column number) to 0, returning true if the element had been 1, false otherwise. 
    inline bool unset(unsigned int,unsigned int);
    /// This method inverts the indicated matrix element (row and column number), i.e. applies the "not" operator to it so 0 -> 1 and 1 -> 0.
    inline void invert(unsigned int,unsigned int);
    /// The overloaded ostream operator pretty prints the matrix, one line per row beginning and ending with square brackets, showing all elements as 0 or 1.
    friend std::ostream& operator <<(std::ostream&,const Binary_Matrix&);
    /// The overloaded unary "not" operator applies this same operator to all the elements of the matrix. 
    friend Binary_Matrix operator !(const Binary_Matrix&);
    /// The overloaded multiplication operator first checks that the two matrices conform to one another (i.e. the number of columns of the first argument is equal to the number of rows of the second), it then carries out the matrix multiplication over GF(2).
    friend Binary_Matrix operator *(const Binary_Matrix&,const Binary_Matrix&);
    /// The overloaded addition operator first checks that the two matrices have the same dimensions and then performs an element-wise "xor" on them to form the sum over GF(2). 
    friend Binary_Matrix operator +(const Binary_Matrix&,const Binary_Matrix&);
    friend class Graph;
  };

  bool Binary_Matrix::get(unsigned int i,unsigned int j) const
  {
    if (i >= nrow) throw std::invalid_argument("The row number argument is illegal for this binary matrix!");
    if (j >= ncolumn) throw std::invalid_argument("The column number argument is illegal for this binary matrix!");
    assert(i >= 0 && j >= 0);
    bool output = (elements[i].count(j) > 0) ? true : false; 
    return output;
  }

  bool Binary_Matrix::set(unsigned int i,unsigned int j) 
  {
    if (i >= nrow) throw std::invalid_argument("The row number argument is illegal for this binary matrix!");
    if (j >= ncolumn) throw std::invalid_argument("The column number argument is illegal for this binary matrix!");
    if (elements[i].count(j) > 0) return false;
    elements[i].insert(j);
    return true;
  }

  bool Binary_Matrix::unset(unsigned int i,unsigned int j)
  {
    if (i >= nrow) throw std::invalid_argument("The row number argument is illegal for this binary matrix!");
    if (j >= ncolumn) throw std::invalid_argument("The column number argument is illegal for this binary matrix!");
    if (elements[i].count(j) > 0) {
      elements[i].erase(j);
      return true;
    }
    return false; 
  }

  void Binary_Matrix::invert(unsigned int i,unsigned int j)
  {
    if (i >= nrow) throw std::invalid_argument("The row number argument is illegal for this binary matrix!");
    if (j >= ncolumn) throw std::invalid_argument("The column number argument is illegal for this binary matrix!");
    if (elements[i].count(j) > 0) {
      elements[i].erase(j);
    }
    else {
      elements[i].insert(j);
    }
  }
}
#endif

