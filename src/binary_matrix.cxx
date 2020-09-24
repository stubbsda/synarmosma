#include "binary_matrix.h"

using namespace SYNARMOSMA;

Binary_Matrix::Binary_Matrix()
{

}

Binary_Matrix::Binary_Matrix(unsigned int n,bool identity)
{
  if (n == 0) throw std::invalid_argument("The matrix dimension must be greater than zero!");
  initialize(n,n);
  if (identity) {
    unsigned int i;
    for(i=0; i<n; ++i) {
      elements[i].insert(i);
    }
  }
}

Binary_Matrix::Binary_Matrix(unsigned int n,unsigned int m,float percent)
{
  if (n == 0 || m == 0) throw std::invalid_argument("The matrix dimensions must be greater than zero!");
  if (percent < -std::numeric_limits<float>::epsilon()) throw std::invalid_argument("The percentage in the binary matrix constructor must not be zero!");
  initialize(n,m);
  if (percent > std::numeric_limits<float>::epsilon()) {
    unsigned int i,j;
    float alpha;
    Random RND;

    for(i=0; i<nrow; ++i) {
      for(j=0; j<ncolumn; ++j) {
        alpha = float(RND.drandom());
        if (alpha < percent) elements[i].insert(j);
      }
    }
  }
}

Binary_Matrix::Binary_Matrix(const Binary_Matrix& source)
{
  unsigned int i;
  nrow = source.nrow;
  ncolumn = source.ncolumn;
  elements = new std::set<unsigned int>[nrow];
  for(i=0; i<nrow; ++i) {
    elements[i] = source.elements[i];
  } 
}

Binary_Matrix& Binary_Matrix::operator =(const Binary_Matrix& source)
{
  if (this == &source) return *this;

  clear();
  unsigned int i;

  nrow = source.nrow;
  ncolumn = source.ncolumn;
  elements = new std::set<unsigned int>[nrow];
  for(i=0; i<nrow; ++i) {
    elements[i] = source.elements[i];
  } 

  return *this;
}

Binary_Matrix::~Binary_Matrix()
{
  if (nrow > 0) delete[] elements;
}

void Binary_Matrix::initialize(unsigned int n,unsigned int m)
{
  if (n == 0 || m == 0) throw std::invalid_argument("The matrix dimensions must be greater than zero!");
  clear();
  nrow = n;
  ncolumn = m;
  elements = new std::set<unsigned int>[nrow];
}

void Binary_Matrix::clear()
{
  if (nrow > 0) delete[] elements;
  nrow = 0;
  ncolumn = 0;
}

int Binary_Matrix::serialize(std::ofstream& s) const
{
  unsigned int i,j,n;
  int count = 0;
  std::set<unsigned int>::const_iterator it;

  s.write((char*)(&nrow),sizeof(int)); count += sizeof(int);
  s.write((char*)(&ncolumn),sizeof(int)); count += sizeof(int);
  for(i=0; i<nrow; ++i) {
    n = elements[i].size();
    s.write((char*)(&n),sizeof(int)); count += sizeof(int);
    for(it=elements[i].begin(); it!=elements[i].end(); ++it) {
      j = *it;
      s.write((char*)(&j),sizeof(int)); count += sizeof(int);
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
  elements = new std::set<unsigned int>[nrow];
  for(i=0; i<nrow; ++i) {
    s.read((char*)(&n),sizeof(int)); count += sizeof(int);
    for(j=0; j<n; ++j) {
      s.read((char*)(&k),sizeof(int)); count += sizeof(int);
      elements[i].insert(k);
    }
  }
  return count;
}

void Binary_Matrix::get_row(unsigned int n,bool* output) const
{
  if (n >= nrow) throw std::invalid_argument("The row number argument is illegal for this binary matrix!");
  unsigned int i;
  std::set<unsigned int>::const_iterator it;

  for(i=0; i<ncolumn; ++i) {
    output[i] = false;
  }
  for(it=elements[n].begin(); it!=elements[n].end(); ++it) {
    output[*it] = true;
  }
}

void Binary_Matrix::multiply(const std::vector<double>& x,std::vector<double>& output) const
{
  if (nrow != x.size()) throw std::invalid_argument("The vector must have a length equal to the number of binary matrix rows!");
  unsigned int i;
  double sum;
  std::set<unsigned int>::const_iterator it;

  output.clear();
  for(i=0; i<nrow; ++i) {
    sum = 0.0;
    for(it=elements[i].begin(); it!=elements[i].end(); ++it) {
      sum += x[*it];
    }
    output.push_back(sum);
  }
}

bool Binary_Matrix::symmetric() const
{
  // This method tests if the matrix is symmetric...
  if (nrow != ncolumn) return false;

  unsigned int i,j;
  std::set<unsigned int> nzero;
  std::set<unsigned int>::const_iterator it;

  for(i=0; i<ncolumn; ++i) {
    // The set nzero will contain all the elements of the i-th column
    for(j=0; j<nrow; ++j) {
      for(it=elements[j].begin(); it!=elements[i].end(); ++it) {
        if (*it == i) nzero.insert(*it); 
      }
    }
    // Is the i-th column the same as the i-th row?
    if (nzero != elements[i]) return false;
    nzero.clear();
  }
  return true;
}

double Binary_Matrix::density() const
{
  unsigned int i,ncount = 0;
  for(i=0; i<nrow; ++i) {
    ncount += elements[i].size();
  }
  double output = double(ncount)/double(nrow*ncolumn);
  return output;
}

int Binary_Matrix::rank() const
{
  unsigned int i,j,k,l,r = 0;
  unsigned int cvector[ncolumn];
  bool found;
  std::set<unsigned int> holder;
  std::set<unsigned int>::const_iterator it;
  Binary_Matrix wcopy(*this);

  for(i=0; i<ncolumn; ++i) {
    cvector[i] = 0;
  }

  for(k=0; k<ncolumn; ++k) {
    found = false;
    for(j=r; j<nrow; ++j) {
      if (wcopy.elements[j].count(k) > 0) {
        found = true;
        break;
      }
    }
    if (!found) continue;
    // Swap rows "j" and "r"...
    holder = wcopy.elements[j];
    wcopy.elements[j] = wcopy.elements[r];
    wcopy.elements[r] = holder;

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i,l,it) firstprivate(cvector) schedule(dynamic,1)
#endif
    for(i=0; i<nrow; ++i) {
      if (wcopy.elements[i].count(k) == 0 || i == r) continue;
      for(it=wcopy.elements[i].begin(); it!=wcopy.elements[i].end(); ++it) {
        cvector[*it] = 1;
      }
      for(it=wcopy.elements[r].begin(); it!=wcopy.elements[r].end(); ++it) {
        l = *it;
        cvector[l] = (cvector[l] == 1) ? 0 : 1;
      }
      wcopy.elements[i].clear();
      for(l=0; l<ncolumn; ++l) {
        if (cvector[l] == 1) {
          wcopy.elements[i].insert(l);
          cvector[l] = 0;
        }
      }
    }
    r++;
  }

  return (signed) r;
}

namespace SYNARMOSMA {
  std::ostream& operator <<(std::ostream& s,const Binary_Matrix& A)
  {
    unsigned int i,j;

    for(i=0; i<A.nrow; ++i) {
      s << "[ ";
      for(j=0; j<A.ncolumn-1; ++j) {
        if (A.elements[i].count(j) == 0) {
          s << "0, ";
        }
        else {
          s << "1, ";
        }
      }
      if (A.elements[i].count(A.ncolumn-1) == 0) {
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

  Binary_Matrix operator !(const Binary_Matrix& A)
  {
    unsigned int i,j;
    Binary_Matrix output(A.nrow,A.ncolumn);

    for(i=0; i<A.nrow; ++i) {
      for(j=0; j<A.ncolumn; ++j) {
        if (A.elements[i].count(j) == 0) output.elements[i].insert(j);
      }
    }
    return output;
  }

  Binary_Matrix operator +(const Binary_Matrix& A,const Binary_Matrix& B)
  {
    if (A.nrow != B.nrow || A.ncolumn != B.ncolumn) throw std::invalid_argument("These two binary matrices do not conform for addition!");

    unsigned int i,j;
    std::set<unsigned int>::const_iterator it;
    Binary_Matrix output(A.nrow,A.ncolumn);

    output = A;

    // Now the matrix addition itself...
    for(i=0; i<B.nrow; ++i) {
      for(it=B.elements[i].begin(); it!=B.elements[i].end(); ++it) {
        j = *it;
        if (output.elements[i].count(j) == 0) {
          output.elements[i].insert(j);
        }
        else {
          // 1 + 1 = 0 over GF(2), so delete this element
          output.elements[i].erase(j);
        }
      }
    }
    return output;
  }

  Binary_Matrix operator *(const Binary_Matrix& A,const Binary_Matrix& B)
  {
    if (A.ncolumn != B.nrow) throw std::invalid_argument("These two binary matrices do not conform for multiplication!");

    unsigned int i,j,l,in1,sum;
    std::set<unsigned int>::const_iterator it;
    Binary_Matrix output(A.nrow,B.ncolumn);

    // Convert the matrix B to a column-oriented format...
    std::vector<unsigned int>* nonzero_row = new std::vector<unsigned int>[B.ncolumn];
    for(l=0; l<B.ncolumn; ++l) {
      for(i=0; i<B.nrow; ++i) {
        if (B.elements[i].count(l) > 0) nonzero_row[l].push_back(i);
      }
    }
    // Now the matrix multiplication itself...
    for(i=0; i<A.nrow; ++i) {
      for(j=0; j<B.ncolumn; ++j) {
        sum = 0;
        for(it=A.elements[i].begin(); it!=A.elements[i].end(); ++it) {
          in1 = *it;
          for(l=0; l<nonzero_row[j].size(); ++l) {
            if (nonzero_row[j][l] == in1) sum += 1;
          }
        }
        if (sum % 2 == 1) output.elements[i].insert(j);
      }
    }
    delete[] nonzero_row;
    return output;
  }
}

