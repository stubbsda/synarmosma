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

#include "binary_matrix.h"

using namespace SYNARMOSMA;

Binary_Matrix::Binary_Matrix()
{
  nrow = 10;
  ncolumn = 10;
  elements = new std::vector<unsigned int>[nrow];
}

Binary_Matrix::Binary_Matrix(unsigned int n)
{
  nrow = n;
  ncolumn = n;
  elements = new std::vector<unsigned int>[nrow];
}

Binary_Matrix::Binary_Matrix(unsigned int n,unsigned int m)
{
  nrow = n;
  ncolumn = m;
  elements = new std::vector<unsigned int>[nrow];
}

Binary_Matrix::Binary_Matrix(const Binary_Matrix& source)
{
  unsigned int i;
  nrow = source.nrow;
  ncolumn = source.ncolumn;
  elements = new std::vector<unsigned int>[nrow];
  for(i=0; i<nrow; ++i) {
    elements[i] = source.elements[i];
  } 
}

Binary_Matrix& Binary_Matrix::operator =(const Binary_Matrix& source)
{
  if (this == &source) return *this;

  unsigned int i;
  delete[] elements;
  nrow = source.nrow;
  ncolumn = source.ncolumn;
  elements = new std::vector<unsigned int>[nrow];
  for(i=0; i<nrow; ++i) {
    elements[i] = source.elements[i];
  } 

  return *this;
}

Binary_Matrix::~Binary_Matrix()
{
  delete[] elements;
}

void Binary_Matrix::initialize(unsigned int n,unsigned int m)
{
  delete[] elements;
  nrow = n;
  ncolumn = m;
  elements = new std::vector<unsigned int>[nrow];
}

bool Binary_Matrix::get(unsigned int i,unsigned int j) const
{
  std::vector<unsigned int>::const_iterator it = std::find(elements[i].begin(),elements[i].end(),j);
  bool output = (it == elements[i].end()) ? false : true; 
  return output;
}

void Binary_Matrix::set(unsigned int i,unsigned int j) 
{
  std::vector<unsigned int>::const_iterator it = std::find(elements[i].begin(),elements[i].end(),j);
  if (it == elements[i].end()) elements[i].push_back(j);
}

void Binary_Matrix::unset(unsigned int i,unsigned int j)
{
  std::vector<unsigned int>::const_iterator it = std::find(elements[i].begin(),elements[i].end(),j);
  if (it != elements[i].end()) elements[i].erase(elements[i].begin() + *it); 
}

unsigned int Binary_Matrix::rank() const
{
  unsigned int i,j,k,r = 0;
  bool found;
  std::vector<unsigned int> nzero;
  std::vector<unsigned int>::const_iterator it,jt;
  Binary_Matrix wcopy(nrow,ncolumn);

  // Copy over the matrix into wcopy
  for(i=0; i<nrow; ++i) {
    wcopy.elements[i] = elements[i];
  }

  for(k=0; k<ncolumn; ++k) {
    found = false;
    for(j=r; j<nrow; ++j) {
      it = std::find(wcopy.elements[j].begin(),wcopy.elements[j].end(),k);
      if (it != wcopy.elements[j].end()) {
        found = true;
        break;
      }
    }
    if (found) {
      nzero = wcopy.elements[j];
      wcopy.elements[j] = wcopy.elements[r];
      wcopy.elements[r] = nzero;
      for(i=0; i<r; ++i) {
        jt = std::find(wcopy.elements[i].begin(),wcopy.elements[i].end(),k);
        if (jt == wcopy.elements[i].end()) continue;
        nzero = wcopy.elements[i];
        for(it=wcopy.elements[r].begin(); it!=wcopy.elements[r].end(); ++it) {
          nzero.push_back(*it);
        }
        // Find the columns that appear twice in nzero and delete both...
        wcopy.elements[i].clear();
        for(it=nzero.begin(); it!=nzero.end(); ++it) {
          j = std::count(nzero.begin(),nzero.end(),*it);
          if (j == 2) continue;
          wcopy.elements[i].push_back(*it);
        }
      }

      for(i=r+1; i<nrow; ++i) {
        it = std::find(wcopy.elements[i].begin(),wcopy.elements[i].end(),k);
        if (it == wcopy.elements[i].end()) continue;
        nzero = wcopy.elements[i];
        for(it=wcopy.elements[r].begin(); it!=wcopy.elements[r].end(); ++it) {
          nzero.push_back(*it);
        }
        // Find the columns that appear twice in nzero and delete both...
        wcopy.elements[i].clear();
        for(it=nzero.begin(); it!=nzero.end(); ++it) {
          j = std::count(nzero.begin(),nzero.end(),*it);
          if (j == 2) continue;
          wcopy.elements[i].push_back(*it);
        }
      }
      r++;      
    }
  }
  return r;
}

namespace SYNARMOSMA {
  std::ostream& operator <<(std::ostream& s,const Binary_Matrix& A)
  {
    unsigned int i,j;
    std::vector<unsigned int>::const_iterator it;

    for(i=0; i<A.nrow; ++i) {
      s << "[ ";
      for(j=0; j<A.ncolumn-1; ++j) {
        it = std::find(A.elements[i].begin(),A.elements[i].end(),j);
        if (it == A.elements[i].end()) {
          s << "0, ";
        }
        else {
          s << "1, ";
        }
      }
      it = std::find(A.elements[i].begin(),A.elements[i].end(),A.ncolumn-1);
      if (it == A.elements[i].end()) {
        s << "0 ";
      }
      else {
        s << "1 ";
      }
      if (i == A.nrow-1) {
        s << "]";
      }
      else {
        s << "]" << std::endl;
      }
    }
    return s;
  }
}

