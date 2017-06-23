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
  class Matrix {
   protected:
    std::vector<std::pair<kind,unsigned int> >* elements;
    unsigned int nrow,ncolumn;
    bool normalized;

    void display(std::ostream&) const;
    bool divisible(unsigned int,unsigned int*,unsigned int*,kind*) const;
    void transpose(const Matrix<kind>&);
    void write_elements(std::ofstream&) const;
    void read_elements(std::ifstream&);
   public:
    static const kind zero;
    static const kind neg1;
    static const kind unity;

    Matrix();
    Matrix(unsigned int);
    Matrix(unsigned int,unsigned int);
    Matrix(const Matrix<kind>&);
    ~Matrix();
    void clear();
    void clear(bool);
    void serialize(std::ofstream&) const;
    void deserialize(std::ifstream&);
    void initialize(unsigned int,unsigned int);
    void set(unsigned int,unsigned int,kind,bool = false);
    kind get(unsigned int,unsigned int) const;
    void increment(unsigned int,unsigned int,kind);
    void convert(kind*,char) const;
    bool empty_row(unsigned int) const;
    inline double sparsity() const {return 1.0 - double(number_nonzero())/double(nrow*ncolumn);};
    unsigned int number_nonzero() const; 
    double dispersion() const;
    bool diagonally_dominant() const;
    bool diagonally_dominant(unsigned int) const;
    bool diagonally_dominant(unsigned int,unsigned int) const;
    unsigned int optimize_dominance(std::vector<unsigned int>&);
    kind determinant() const;
    void multiply(const std::vector<kind>&,std::vector<kind>&) const;
    void increment(const Matrix<kind>&);
    void get_diagonal(std::vector<kind>&) const;
    int gauss_seidel_solver(std::vector<kind>&,const std::vector<kind>&,double,int);
    kind get_first_nonzero(unsigned int) const;
    inline int get_nrow() const {return nrow;};
    void get_row(std::vector<kind>&,std::vector<unsigned int>&,int) const;
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
}
#endif
