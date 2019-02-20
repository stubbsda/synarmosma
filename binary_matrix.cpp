#include "binary_matrix.h"

using namespace SYNARMOSMA;

Binary_Matrix::Binary_Matrix()
{
  nrow = 10;
  ncolumn = 10;
  elements = new std::vector<unsigned int>[nrow];
}

Binary_Matrix::Binary_Matrix(int n)
{
  assert(n > 0);
  nrow = n;
  ncolumn = n;
  elements = new std::vector<unsigned int>[nrow];
}

Binary_Matrix::Binary_Matrix(int n,int m)
{
  assert(n > 0 && m > 0);
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

void Binary_Matrix::initialize(int n,int m)
{
  assert(n > 0 && m > 0);
  delete[] elements;
  nrow = n;
  ncolumn = m;
  elements = new std::vector<unsigned int>[nrow];
}

void Binary_Matrix::clear()
{
  nrow = 0;
  ncolumn = 0;
  delete[] elements;
}

int Binary_Matrix::serialize(std::ofstream& s) const
{
  unsigned int i,j,n;
  int count = 0;

  s.write((char*)(&nrow),sizeof(int)); count += sizeof(int);
  s.write((char*)(&ncolumn),sizeof(int)); count += sizeof(int);
  for(i=0; i<nrow; ++i) {
    n = elements[i].size();
    s.write((char*)(&n),sizeof(int)); count += sizeof(int);
    for(j=0; j<n; ++j) {
      s.write((char*)(&elements[i][j]),sizeof(int)); count += sizeof(int);
    }
  }
  return count;
}

int Binary_Matrix::deserialize(std::ifstream& s)
{
  unsigned int i,j,k,n;
  int count = 0;

  clear();

  s.read((char*)(&nrow),sizeof(int)); count += sizeof(int);
  s.read((char*)(&ncolumn),sizeof(int)); count += sizeof(int);
  elements = new std::vector<unsigned int>[nrow];
  for(i=0; i<nrow; ++i) {
    s.read((char*)(&n),sizeof(int)); count += sizeof(int);
    for(j=0; j<n; ++j) {
      s.read((char*)(&k),sizeof(int)); count += sizeof(int);
      elements[i].push_back(k);
    }
  }
  return count;
}

void Binary_Matrix::get_row(int n,bool* output) const
{
  assert(n >= 0);
  unsigned int i;
  std::vector<unsigned int>::const_iterator it;

  for(i=0; i<ncolumn; ++i) {
    output[i] = false;
  }
  for(it=elements[n].begin(); it!=elements[n].end(); ++it) {
    output[*it] = true;
  }
}

bool Binary_Matrix::get(int i,int j) const
{
  assert(i >= 0 && i >= 0);
  std::vector<unsigned int>::const_iterator it = std::find(elements[i].begin(),elements[i].end(),j);
  bool output = (it == elements[i].end()) ? false : true; 
  return output;
}

void Binary_Matrix::set(int i,int j) 
{
  assert(i >= 0 && i >= 0);
  std::vector<unsigned int>::const_iterator it = std::find(elements[i].begin(),elements[i].end(),j);
  if (it == elements[i].end()) elements[i].push_back(j);
}

void Binary_Matrix::unset(int i,int j)
{
  assert(i >= 0 && i >= 0);
  std::vector<unsigned int>::const_iterator it = std::find(elements[i].begin(),elements[i].end(),j);
  if (it != elements[i].end()) elements[i].erase(elements[i].begin() + *it); 
}

int Binary_Matrix::rank() const
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
  return (signed) r;
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

  Binary_Matrix operator +(const Binary_Matrix& A,const Binary_Matrix& B)
  {
#ifdef DEBUG
    assert(A.nrow == B.nrow && A.ncolumn == B.ncolumn);
#endif

    unsigned int i,j,k;
    std::vector<unsigned int>::const_iterator it;
    Binary_Matrix output(A.nrow,A.ncolumn);

    output = A;

    // Now the matrix addition itself...
    for(i=0; i<B.nrow; ++i) {
      for(j=0; j<B.elements[i].size(); ++j) {
        k = B.elements[i][j];
        it = std::find(output.elements[i].begin(),output.elements[i].end(),k);
        if (it == output.elements[i].end()) {
          output.elements[i].push_back(k);
        }
        else {
          // 1 + 1 = 0 over GF(2), so delete this element
          output.elements[i].erase(output.elements[i].begin() + *it);
        }
      }
    }
    return output;
  }

  Binary_Matrix operator *(const Binary_Matrix& A,const Binary_Matrix& B)
  {
#ifdef DEBUG
    assert(A.ncolumn == B.nrow);
#endif
    unsigned int i,j,k,l,in1,sum;
    Binary_Matrix output(A.nrow,B.ncolumn);

    // Convert the matrix B to a column-oriented format...
    std::vector<unsigned int>* nonzero_row = new std::vector<unsigned int>[B.ncolumn];
    for(l=0; l<B.ncolumn; ++l) {
      for(i=0; i<B.nrow; ++i) {
        for(j=0; j<B.elements[i].size(); ++j) {
          if (B.elements[i][j] == l) nonzero_row[l].push_back(i);
        }
      }
    }
    // Now the matrix multiplication itself...
    for(i=0; i<A.nrow; ++i) {
      for(j=0; j<B.ncolumn; ++j) {
        sum = 0;
        for(k=0; k<A.elements[i].size(); ++k) {
          in1 = A.elements[i][k];
          for(l=0; l<nonzero_row[j].size(); ++l) {
            if (nonzero_row[j][l] == in1) sum += 1;
          }
        }
        if (sum % 2 == 1) output.elements[i].push_back(j);
      }
    }
    delete[] nonzero_row;
    return output;
  }
}

