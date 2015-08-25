/*
  Copyright 2014 Daniel Stubbs

  This file is part of Synarmosma.

  Synarmosma is free software: you can redistribute it and/or modify 
  it under the terms of the GNU General Public License as published by 
  the Free Software Foundation, either version 3 of the License, or 
  (at your option) any later version.

  Synarmosma is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with Synarmosma.  If not, see <http://www.gnu.org/licenses/>.
*/

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
   private:
    std::vector<std::pair<kind,unsigned int> >* elements;
    unsigned int nrow,ncolumn;
    bool normalized;

    void display(std::ostream&) const;
    bool divisible(unsigned int,unsigned int*,unsigned int*,kind*) const;
    void transpose(const Matrix<kind>&);
    void clear();
    void clear(bool);
   public:
    static const kind zero;
    static const kind neg1;
    static const kind unity;

    Matrix();
    Matrix(unsigned int);
    Matrix(unsigned int,unsigned int);
    Matrix(const Matrix<kind>&);
    ~Matrix();
    void serialize(std::ofstream&) const;
    void deserialize(std::ifstream&);
    void initialize(unsigned int,unsigned int);
    void set(unsigned int,unsigned int,kind);
    kind get(unsigned int,unsigned int) const;
    bool empty_row(unsigned int) const;
    kind get_first_nonzero(unsigned int) const;
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
