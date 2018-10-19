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
  void permute(Matrix<kind>&,int,int,char);

  template<class kind>
  void invert(Matrix<kind>&,int,char);

  template<class kind>
  void combine(Matrix<kind>&,int,int,const kind&,char);

  template<class kind>
  void pivot_row(Matrix<kind>&,int,int);

  template<class kind>
  void pivot_column(Matrix<kind>&,int,int);

  template<class kind>
  void move_minimum(Matrix<kind>&,int);

  template<class kind>
  void prepare_matrix(Matrix<kind>&,int);

  template<class kind>
  unsigned int normalize(Matrix<kind>&);

  template<class kind>
  class Matrix {
   protected:
    std::vector<std::pair<kind,unsigned int> >* elements;
    unsigned int nrow,ncolumn;
    bool normalized;

    void display(std::ostream&) const;
    bool divisible(int,unsigned int*,unsigned int*,kind*) const;
    void transpose(const Matrix<kind>&);
    int write_elements(std::ofstream&) const;
    int read_elements(std::ifstream&);
   public:
    static const kind zero;
    static const kind neg1;
    static const kind unity;

    Matrix();
    Matrix(int);
    Matrix(int,int);
    Matrix(const Matrix<kind>&);
    ~Matrix();
    void clear();
    void clear(bool);
    int serialize(std::ofstream&) const;
    int deserialize(std::ifstream&);
    void initialize(int,int);
    inline void set(int,int,kind,bool = false);
    inline kind get(int,int) const;
    void increment(int,int,kind);
    void convert(kind*,char) const;
    bool empty_row(int) const;
    inline double sparsity() const {return 1.0 - double(number_nonzero())/double(nrow*ncolumn);};
    int number_nonzero() const; 
    double dispersion() const;
    kind diagonal_element(int) const;
    void get_diagonal(std::vector<kind>&) const;
    bool diagonally_dominant() const;
    bool diagonally_dominant(int) const;
    bool diagonally_dominant(int,int) const;
    int optimize_dominance(std::vector<unsigned int>&);
    kind determinant() const;
    void multiply(const std::vector<kind>&,std::vector<kind>&) const;
    void increment(const Matrix<kind>&);
    void homotopy_scaling(kind,kind,Matrix<kind>*) const;
    void gauss_seidel_solver(std::vector<kind>&,const std::vector<kind>&,double,int);
    kind get_first_nonzero(int) const;
    inline int get_nrow() const {return nrow;};
    void get_row(std::vector<kind>&,std::vector<unsigned int>&,int) const;
    Matrix<kind>& operator =(const Matrix<kind>&);
    friend void permute<>(Matrix<kind>&,int,int,char);
    friend void invert<>(Matrix<kind>&,int,char);
    friend void combine<>(Matrix<kind>&,int,int,const kind&,char);
    friend void pivot_row<>(Matrix<kind>&,int,int);
    friend void pivot_column<>(Matrix<kind>&,int,int);
    friend void move_minimum<>(Matrix<kind>&,int);
    friend void prepare_matrix<>(Matrix<kind>&,int);
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
  kind Matrix<kind>::get(int n,int m) const
  {
    assert(n >= 0 && m >= 0);
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
  void Matrix<kind>::set(int n,int m,kind v,bool increment)
  {
    assert(n >= 0 && m >= 0);
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
}
#endif
