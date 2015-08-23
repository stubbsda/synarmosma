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

#ifndef _bmatrixh
#define _bmatrixh

namespace SYNARMOSMA {
  class Binary_Matrix {
   private:
    std::vector<unsigned int>* elements;
    unsigned int nrow,ncolumn;

   public:
    Binary_Matrix();
    Binary_Matrix(unsigned int);
    Binary_Matrix(unsigned int,unsigned int);
    Binary_Matrix(const Binary_Matrix&);
    Binary_Matrix& operator =(const Binary_Matrix&);
    ~Binary_Matrix();
    void initialize(unsigned int,unsigned int);
    void get_row(unsigned int,bool*) const;
    bool get(unsigned int,unsigned int) const;
    void set(unsigned int,unsigned int);
    void unset(unsigned int,unsigned int);
    unsigned int rank() const;
    friend std::ostream& operator <<(std::ostream&,const Binary_Matrix&);
    friend Binary_Matrix operator *(const Binary_Matrix&,const Binary_Matrix&);
    friend Binary_Matrix operator +(const Binary_Matrix&,const Binary_Matrix&);
    friend class Graph;
  };
}
#endif

