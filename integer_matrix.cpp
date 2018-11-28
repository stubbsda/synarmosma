#include "integer_matrix.h"

using namespace SYNARMOSMA;

extern Random RND;

template<class kind>
Integer_Matrix<kind>::Integer_Matrix()
{

}

template<class kind>
Integer_Matrix<kind>::Integer_Matrix(unsigned int n,bool identity)
{
  if (n == 0) throw std::invalid_argument("The matrix dimension must be greater than zero!");
  nrow = n;
  ncolumn = n;
  elements = new std::vector<std::pair<kind,unsigned int> >[nrow];
  if (identity) { 
    unsigned int i;
    for(i=0; i<nrow; ++i) {
      elements[i].push_back(std::pair<kind,unsigned int>(Integer_Matrix<kind>::unity,i));
    }
  }
}

template<class kind>
Integer_Matrix<kind>::Integer_Matrix(unsigned int n,unsigned int m)
{
  if (n == 0 || m == 0) throw std::invalid_argument("The matrix dimensions must be greater than zero!");
  nrow = n;
  ncolumn = m;
  elements = new std::vector<std::pair<kind,unsigned int> >[nrow];
}

template<class kind>
Integer_Matrix<kind>::Integer_Matrix(const Integer_Matrix<kind>& source)
{
  unsigned int i;

  clear();
  nrow = source.nrow;
  ncolumn = source.ncolumn;
  normalized = source.normalized;
  elements = new std::vector<std::pair<kind,unsigned int> >[nrow];
  for(i=0; i<nrow; ++i) {
    elements[i] = source.elements[i];
  }
}

template<class kind>
Integer_Matrix<kind>& Integer_Matrix<kind>::operator =(const Integer_Matrix<kind>& source)
{
  if (this == &source) return *this;

  unsigned int i;

  clear();
  nrow = source.nrow;
  ncolumn = source.ncolumn;
  normalized = source.normalized;
  elements = new std::vector<std::pair<kind,unsigned int> >[nrow];
  for(i=0; i<nrow; ++i) {
    elements[i] = source.elements[i];
  }

  return *this;
}

template<class kind>
Integer_Matrix<kind>::~Integer_Matrix()
{
  if (nrow > 0) delete[] elements;
}

template<class kind>
void Integer_Matrix<kind>::display(std::ostream& s) const
{
  unsigned int i,j,k;
  kind w;
  bool found;

  for(i=0; i<nrow; ++i) {
    for(j=0; j<ncolumn; ++j) {
      w = Integer_Matrix<kind>::zero;
      found = false;
      for(k=0; k<elements[i].size(); ++k) {
        if (elements[i][k].second == j) {
          w = elements[i][k].first;
          found = true;
          break;
        }
      }
      if (!found) {
        s << "%" << w << " ";
      }
      else {
        s << w << " ";
      }
    }
    s << std::endl;
  }
  s << std::endl;
}

template<class kind>
void Integer_Matrix<kind>::convert(kind* A,char mtype) const
{
  if (mtype != 'r' && mtype != 'c') throw std::invalid_argument("The conversion type for the matrix must be row-oriented or column-oriented!");
  // This method supposes that A has already been allocated the appropriate amount
  // of memory...
  unsigned int i,j;

  for(i=0; i<nrow; ++i) {
    for(j=0; j<ncolumn; ++j) {
      A[nrow*i+j] = Integer_Matrix<kind>::zero;
    }
  }
  if (mtype == 'r') {
    // Row major format
    for(i=0; i<nrow; ++i) {
      for(j=0; j<elements[i].size(); ++j) {
        A[nrow*i+elements[i][j].second] = elements[i][j].first;
      }
    }
  }
  else {
    // Column major format
    for(i=0; i<nrow; ++i) {
      for(j=0; j<elements[i].size(); ++j) {
        A[ncolumn*elements[i][j].second + i] = elements[i][j].first;
      }
    }
  }  
}

template<class kind>
kind Integer_Matrix<kind>::get(unsigned int n,unsigned int m) const
{
  if (n >= nrow) throw std::invalid_argument("The row number argument is illegal for this matrix!");
  if (m >= ncolumn) throw std::invalid_argument("The column number argument is illegal for this matrix!");
  unsigned int i,mu = m;
  kind output = Integer_Matrix<kind>::zero; 

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
kind Integer_Matrix<kind>::get_diagonal(unsigned int n) const
{
  if (n >= nrow) throw std::invalid_argument("The row number argument is illegal for this matrix!");
  unsigned int i;
  kind output = Integer_Matrix<kind>::zero;

  for(i=0; i<elements[n].size(); ++i) {
    if (elements[n][i].second == n) {
      output = elements[n][i].first;
      break;
    }
  }
  return output;
}

template<class kind>
bool Integer_Matrix<kind>::divisible(unsigned int n,unsigned int* out1,unsigned int* out2,kind* f) const
{
  if (n >= nrow) throw std::invalid_argument("The row number argument is illegal for this matrix!");
  unsigned int i,j;
  kind b1,b2 = Integer_Matrix<kind>::zero;

  for(i=0; i<elements[n].size(); ++i) {
    if (elements[n][i].second == n) b2 = elements[n][i].first;
  }
  if (b2 == Integer_Matrix<kind>::zero) return false;
  for(i=n+1; i<nrow; ++i) {
    for(j=0; j<elements[i].size(); ++j) {
      if (elements[i][j].second <= n) continue;
      b1 = elements[i][j].first;
      if ((b1 % b2) != 0) {
        *out1 = i;
        *out2 = elements[i][j].second;
        *f = b1/b2;
        return true;
      }
    }
  }
  return false;
}

template<class kind>
void Integer_Matrix<kind>::increment(const Integer_Matrix<kind>& arg)
{
  if (nrow != arg.nrow || ncolumn != arg.ncolumn) throw std::invalid_argument("The dimensions of the two matrices must be identical to add them!");
  unsigned int i,j,k,n;
  bool found;

  for(i=0; i<nrow; ++i) {
    for(j=0; j<arg.elements[i].size(); ++j) {
      found = false;
      n = arg.elements[i][j].second;
      for(k=0; k<elements[i].size(); ++k) {
        if (elements[i][k].second == n) {
          elements[i][k].first += arg.elements[i][j].first;
          found = true;
          break;
        }
      }
      if (!found) elements[i].push_back(arg.elements[i][j]); 
    }
  }
}

template<class kind>
void Integer_Matrix<kind>::multiply(const std::vector<kind>& b,std::vector<kind>& output) const
{
  if (ncolumn != b.size()) throw std::invalid_argument("The vector's length must equal the number of matrix columns!");
  unsigned int i,j;
  kind sum;

  output.clear();

  for(i=0; i<nrow; ++i) {
    sum = Integer_Matrix<kind>::zero;
    for(j=0; j<elements[i].size(); ++j) {
      sum += elements[i][j].first*b[elements[i][j].second];
    }
    output.push_back(sum);
  }
}

template<class kind>
kind Integer_Matrix<kind>::determinant() const
{
  if (nrow != ncolumn) throw std::runtime_error("The determinant can only be computed for square matrices!");
  unsigned int i,j,k,q;
  bool found;
  std::vector<int> symbols,rest;
  std::vector<std::vector<int> > plist;
  kind product,output = Integer_Matrix<kind>::zero;

  for(i=0; i<nrow; ++i) {
    rest.push_back(i);
  }
  permute(plist,symbols,rest,nrow);

  for(i=0; i<plist.size(); ++i) {
    product = (parity(plist[i]) == 1) ? Integer_Matrix<kind>::unity : Integer_Matrix<kind>::neg1;
    for(j=0; j<nrow; ++j) {
      q = (unsigned) plist[i][j];
      // Does the j-th row contain an element in column q?
      found = false;
      for(k=0; k<elements[j].size(); ++k) {
        if (elements[j][k].second == q) {
          product *= elements[j][k].first;
          found = true;
          break;
        }
      }
      if (!found) {
        product = Integer_Matrix<kind>::zero;
        break;
      }
    }
    output += product;
  }

  return output;
}

template<class kind>
void Integer_Matrix<kind>::get_diagonal(std::vector<kind>& output) const
{
  unsigned int i,j;
  kind d;

  output.clear();

  for(i=0; i<nrow; ++i) {
    d = Integer_Matrix<kind>::zero;
    for(j=0; j<elements[i].size(); ++j) {
      if (elements[i][j].second == i) {
        d = elements[i][j].first;
        break;
      }
    }
    output.push_back(d);
  }
}

template<class kind>
void Integer_Matrix<kind>::clear(bool deep)
{
  normalized = false;
  if (deep) {
    if (nrow > 0) delete[] elements;
    nrow = 0;
    ncolumn = 0;
  }
  else {
    unsigned int i;
    for(i=0; i<nrow; ++i) {
      elements[i].clear();
    }
  }
}

template<class kind>
bool Integer_Matrix<kind>::consistent() const
{
  if (nrow == 0 || ncolumn == 0) return false;
  unsigned int i,j,m;
  std::set<unsigned int> indices;

  for(i=0; i<nrow; ++i) {
    for(j=0; j<elements[i].size(); ++j) {
      m = elements[i][j].second;
      // Make sure the column index lies in the appropriate range...
      if (m >= ncolumn) return false;
      // Also make sure that the there is only element corresponding to this row and column...
      if (indices.count(m) > 0) return false;
      indices.insert(m);
    }
    indices.clear();
  }
  return true;
}

namespace SYNARMOSMA {
  template<>
  int Integer_Matrix<NTL::ZZ>::write_elements(std::ofstream& s) const
  {
    unsigned int i,j,k,n;
    int count = 0;

    for(i=0; i<nrow; ++i) {
      n = elements[i].size();
      s.write((char*)(&n),sizeof(int)); count += sizeof(int);
      for(j=0; j<n; ++j) {        
        count += write_ZZ(s,elements[i][j].first);
        k = elements[i][j].second;
        s.write((char*)(&k),sizeof(int)); count += sizeof(int);
      }
    }
    return count;
  }
}

template<class kind>
int Integer_Matrix<kind>::write_elements(std::ofstream& s) const
{
  unsigned int i,j,k,n;
  int count = 0;
  kind x;

  for(i=0; i<nrow; ++i) {
    n = elements[i].size();
    s.write((char*)(&n),sizeof(int)); count += sizeof(int);
    for(j=0; j<n; ++j) {
      x = elements[i][j].first;
      k = elements[i][j].second;
      s.write((char*)(&x),sizeof(kind)); count += sizeof(kind);
      s.write((char*)(&k),sizeof(int)); count += sizeof(int);
    }
  }
  return count;
}

template<class kind>
int Integer_Matrix<kind>::serialize(std::ofstream& s) const
{
  int count = 0;

  s.write((char*)(&nrow),sizeof(int)); count += sizeof(int);
  s.write((char*)(&ncolumn),sizeof(int)); count += sizeof(int);
  s.write((char*)(&normalized),sizeof(bool)); count += sizeof(bool);
  count += write_elements(s);

  return count;
}

namespace SYNARMOSMA {
  template<>
  int Integer_Matrix<NTL::ZZ>::read_elements(std::ifstream& s)
  {
    unsigned int i,j,k,n;
    int count = 0;
    NTL::ZZ p;

    for(i=0; i<nrow; ++i) {
      s.read((char*)(&n),sizeof(int)); count += sizeof(int);
      for(j=0; j<n; ++j) {
        count += read_ZZ(s,p);
        s.read((char*)(&k),sizeof(int)); count += sizeof(int);
        elements[i].push_back(std::pair<NTL::ZZ,unsigned int>(p,k));
      }
    }
    return count;
  }
}

template<class kind>
int Integer_Matrix<kind>::read_elements(std::ifstream& s)
{
  unsigned int i,j,k,n;
  int count = 0;
  kind x;

  for(i=0; i<nrow; ++i) {
    s.read((char*)(&n),sizeof(int)); count += sizeof(int);
    for(j=0; j<n; ++j) {
      s.read((char*)(&x),sizeof(kind)); count += sizeof(kind);
      s.read((char*)(&k),sizeof(int)); count += sizeof(int);
      elements[i].push_back(std::pair<kind,unsigned int>(x,k));
    }
  }
  return count;
}

template<class kind>
int Integer_Matrix<kind>::deserialize(std::ifstream& s)
{
  int count = 0;

  clear();

  s.read((char*)(&nrow),sizeof(int)); count += sizeof(int);
  s.read((char*)(&ncolumn),sizeof(int)); count += sizeof(int);
  s.read((char*)(&normalized),sizeof(bool)); count += sizeof(bool);
  elements = new std::vector<std::pair<kind,unsigned int> >[nrow];
  count += read_elements(s);

  return count;
}

template<class kind>
unsigned int Integer_Matrix<kind>::eliminate_zeros()
{
  // This method finds ane eliminates zeros from the matrix elements...
  unsigned int i,j,nmodify = 0;
  std::vector<std::pair<kind,unsigned int> > rvector;

  for(i=0; i<nrow; ++i) {
    rvector.clear();
    for(j=0; j<elements[i].size(); ++j) {
      if (elements[i][j].first == Integer_Matrix<kind>::zero) {
        nmodify++;
        continue;
      }
      rvector.push_back(elements[i][j]);
    }
    if (rvector.size() != elements[i].size()) elements[i] = rvector;
  }
  return nmodify;
}

template<class kind>
kind Integer_Matrix<kind>::get_first_nonzero(unsigned int r) const
{
  if (r >= nrow) throw std::invalid_argument("The row number argument is illegal for this matrix!");
  if (elements[r].empty()) return Integer_Matrix<kind>::zero;
  unsigned int i,min_index = ncolumn;
  kind output = Integer_Matrix<kind>::zero;

  for(i=0; i<elements[r].size(); ++i) {
    if (elements[r][i].second < min_index && elements[r][i].first != Integer_Matrix<kind>::zero) {
      output = elements[r][i].first;
      min_index = elements[r][i].second;
    }
  }
  return output; 
}

template<class kind>
unsigned int Integer_Matrix<kind>::number_nonzero() const
{
  unsigned int i,j,nzero = 0;

  for(i=0; i<nrow; ++i) {
    for(j=0; j<elements[i].size(); ++j) {
      if (elements[i][j].first != Integer_Matrix<kind>::zero) nzero++;
    }
  }
  return nzero;
}

template<class kind>
bool Integer_Matrix<kind>::empty_row(unsigned int r) const
{
  if (r >= nrow) throw std::invalid_argument("The row number argument is illegal for this matrix!");
  unsigned int i;
  bool output = true;

  for(i=0; i<elements[r].size(); ++i) {
    if (elements[r][i].first != Integer_Matrix<kind>::zero) {
      output = false;
      break;
    }
  }
  return output;
}

template<class kind>
void Integer_Matrix<kind>::initialize(unsigned int n,unsigned int m)
{
  if (n == 0 || m == 0) throw std::invalid_argument("The matrix dimensions must be greater than zero!");
  if (nrow == n) {
    clear(false);
  }
  else {
    if (nrow > 0) clear(true);
    nrow = n;
    elements = new std::vector<std::pair<kind,unsigned int> >[nrow];
  }
  ncolumn = m;
}

template<class kind>
void Integer_Matrix<kind>::transpose(const Integer_Matrix<kind>& source)
{
  unsigned int i,j;
  std::pair<kind,unsigned int> p;
  delete[] elements;

  nrow = source.ncolumn;
  ncolumn = source.nrow;
  elements = new std::vector<std::pair<kind,unsigned int> >[nrow];
  for(i=0; i<source.nrow; ++i) {
    for(j=0; j<source.elements[i].size(); ++j) {
      p.first = source.elements[i][j].first;
      p.second = i;
      elements[source.elements[i][j].second].push_back(p);
    }
  }
}

template<class kind>
std::vector<kind>& operator *(const Integer_Matrix<kind>& A,const std::vector<kind>& b)
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
Integer_Matrix<kind>& operator +(const Integer_Matrix<kind>& A,const Integer_Matrix<kind>& B)
{
  if (A.nrow != B.nrow || A.ncolumn != B.ncolumn) throw std::invalid_argument("The dimensions of the two matrices don't conform for addition!");
  unsigned int i,j,k,in1,cvalue;
  bool found;
  Integer_Matrix<kind> output(A.nrow,A.ncolumn);

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
Integer_Matrix<kind>& operator *(const Integer_Matrix<kind>& A,const Integer_Matrix<kind>& B)
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
void Integer_Matrix<kind>::get_row(std::vector<kind>& v,std::vector<unsigned int>& c,unsigned int r) const
{
  if (r >= nrow) throw std::invalid_argument("The row number argument is illegal for this matrix!");
  unsigned int i;
  v.clear();
  c.clear();
  for(i=0; i<elements[r].size(); ++i) {
    v.push_back(elements[r][i].first);
    c.push_back(elements[r][i].second);
  }
}

namespace SYNARMOSMA {
  template<class kind>
  void permute(Integer_Matrix<kind>& A,unsigned int n,unsigned int m,char type)
  {
    // If (type == c), then column(n) <-> column(m)
    // If (type == r), then row(n) <-> row(m)
    if (type != 'c' && type != 'r') throw std::invalid_argument("The matrix permutation type must be row-based or column-based!");
    if (n == m) throw std::invalid_argument("The row or column indices to be permuted must be distinct!");
    if (type == 'r' && (n >= A.nrow || m >= A.nrow)) throw std::invalid_argument("The row number argument is illegal for this matrix!");
    if (type == 'c' && (n >= A.ncolumn || m >= A.ncolumn)) throw std::invalid_argument("The column number argument is illegal for this matrix!");

    if (type == 'c') {
      unsigned int i,j;
      for(i=0; i<A.nrow; ++i) {
        for(j=0; j<A.elements[i].size(); ++j) {
          if (A.elements[i][j].second == n) {
            A.elements[i][j].second = m;
          }
          else if (A.elements[i][j].second == m) {
            A.elements[i][j].second = n;
          }
        }
      }
    }
    else {
      std::vector<std::pair<kind,unsigned int> > pr = A.elements[n];
      A.elements[n] = A.elements[m];
      A.elements[m] = pr;
    }
  }

  template<class kind>
  void invert(Integer_Matrix<kind>& A,unsigned int n,char type)
  {
    // If (type == c), then column(n) -> -1*column(n)
    // If (type == r), then row(n) -> -1*row(n)
    if (type != 'c' && type != 'r') throw std::invalid_argument("The matrix inversion type must be row-based or column-based!");
    if (type == 'r' && n >= A.nrow) throw std::invalid_argument("The row number argument is illegal for this matrix inversion!");
    if (type == 'c' && n >= A.ncolumn) throw std::invalid_argument("The column number argument is illegal for this matrix inversion!");
    unsigned int i,j;

    if (type == 'c') {
      for(i=0; i<A.nrow; ++i) {
        for(j=0; j<A.elements[i].size(); ++j) {
          if (A.elements[i][j].second == n) A.elements[i][j].first *= Integer_Matrix<kind>::neg1;
        }
      }
    }
    else {
      for(i=0; i<A.elements[n].size(); ++i) {
        A.elements[n][i].first *= Integer_Matrix<kind>::neg1;
      }
    }
  }

  template<class kind>
  void combine(Integer_Matrix<kind>& A,unsigned int n,unsigned int m,const kind& q,char type)
  {
    // If (type == r), then row(n) -> row(n) + q*row(m)
    // If (type == c), then column(n) -> column(n) + q*column(m)
    if (type != 'c' && type != 'r') throw std::invalid_argument("The matrix combination type must be row-based or column-based!");
    if (type == 'r' && (n >= A.nrow || m >= A.nrow)) throw std::invalid_argument("The row number argument is illegal for this matrix!");
    if (type == 'c' && (n >= A.ncolumn || m >= A.ncolumn)) throw std::invalid_argument("The column number argument is illegal for this matrix!");
    unsigned int i,j,l = 0;
    bool found,fd1,fd2;
    kind wv = Integer_Matrix<kind>::zero;

    if (type == 'c') {
      for(i=0; i<A.nrow; ++i) {
        fd1 = false;
        fd2 = false;
        for(j=0; j<A.elements[i].size(); ++j) {
          if (A.elements[i][j].second == n) {
            fd1 = true;
            wv = q*A.elements[i][j].first;
            if (fd2) break;
          }
          else if (A.elements[i][j].second == m) {
            fd2 = true;
            l = j;
            if (fd1) break;
          }
        }
        if (!fd1) continue;
        if (!fd2) {
          A.elements[i].push_back(std::pair<kind,unsigned int>(wv,m));
        }
        else {
          A.elements[i][l].first += wv;
        }
      }
    }
    else {
      for(i=0; i<A.elements[m].size(); ++i) {
        l = A.elements[m][i].second;
        found = false;
        for(j=0; j<A.elements[n].size(); ++j) {
          if (A.elements[n][j].second == l) {
            A.elements[n][j].first += q*A.elements[m][i].first;
            found = true;
            break;
          }
        }
        if (!found) A.elements[n].push_back(std::pair<kind,unsigned int>(q*A.elements[m][i].first,l));
      }
    }
    A.eliminate_zeros();
  }

  template<class kind>
  void pivot_row(Integer_Matrix<kind>& A,unsigned int k,unsigned int l)
  {
    if (k >= A.nrow) throw std::invalid_argument("The row number argument is illegal for this matrix pivot!");
    if (l >= A.ncolumn) throw std::invalid_argument("The column number argument is illegal for this matrix pivot!");
    unsigned int i,j;
    bool found;
    kind n,q,d = Integer_Matrix<kind>::zero;

    for(i=0; i<A.elements[k].size(); ++i) {
      if (A.elements[k][i].second == l) {
        d = A.elements[k][i].first;
        break;
      }
    }
    for(i=k+1; i<A.nrow; ++i) {
      found = false;
      for(j=0; j<A.elements[i].size(); ++j) {
        if (A.elements[i][j].second == l) {
          n = A.elements[i][j].first;
          found = true;
          break;
        }
      }
      if (!found) continue;
      q = n/d; 
      combine(A,i,k,-q,'r');
    }
  }

  template<class kind>
  void pivot_column(Integer_Matrix<kind>& A,unsigned int k,unsigned int l)
  {
    if (k >= A.nrow) throw std::invalid_argument("The row number argument is illegal for this matrix!");
    if (l >= A.ncolumn) throw std::invalid_argument("The column number argument is illegal for this matrix!");
    unsigned int i,j;
    bool found;
    kind n,q,d = Integer_Matrix<kind>::zero;
 
    for(i=0; i<A.elements[k].size(); ++i) {
      if (A.elements[k][i].second == l) {
        d = A.elements[k][i].first;
        break;
      }
    }
    for(i=l+1; i<A.ncolumn; ++i) {
      found = false;
      for(j=0; j<A.elements[k].size(); ++j) {
        if (A.elements[k][j].second == i) {
          n = A.elements[k][j].first;
          found = true;
          break;
        }
      }
      if (!found) continue;
      q = n/d;
      combine(A,k,i,-q,'c');
    }
  }

  template<class kind>
  void move_minimum(Integer_Matrix<kind>& A,unsigned int n)
  {
    if (n >= A.nrow) throw std::invalid_argument("The row number argument is illegal for this matrix!");
    unsigned int i,j,im,jm,istart = 0,ulimit = A.nrow - n;  
    kind beta,alpha = Integer_Matrix<kind>::zero;
    bool full[ulimit];
    std::pair<kind,unsigned int>* VY = new std::pair<kind,unsigned int>[ulimit];
    std::vector<std::pair<kind,unsigned int> > VX;

    for(i=n; i<A.nrow; ++i) {
      full[i-n] = false;
      if (A.empty_row(i)) continue;
      for(j=0; j<A.elements[i].size(); ++j) {
        if (A.elements[i][j].second >= n) {
          VX.push_back(A.elements[i][j]);
        }
      }
      full[i-n] = true;
      if (VX[0].first < 0) {
        alpha = -VX[0].first;
      }
      else {
        alpha = VX[0].first;
      }
      jm = VX[0].second;
      for(j=1; j<VX.size(); ++j) {
        if (VX[j].first < 0) {
          beta = -VX[j].first;
        }
        else {
          beta = VX[j].first;
        }
        if (beta < alpha) {
          alpha = beta;
          jm = VX[j].second;
        }
      }
      VY[i-n] = std::pair<kind,unsigned int>(alpha,jm);
      VX.clear();
    }
    for(i=0; i<ulimit; ++i) {
      if (!full[i]) continue;
      alpha = VY[i].first;
      istart = i;
      break;
    }
    im = istart;
    for(i=istart+1; i<ulimit; ++i) {
      if (!full[i]) continue;
      beta = VY[i].first;
      if (beta < alpha) {
        alpha = beta;
        im = i;
      }
    }
    jm = VY[im].second;
    im += n;
    if (n != im) permute(A,n,im,'r');
    if (n != jm) permute(A,n,jm,'c');

    delete[] VY;
  }

  template<class kind>
  void prepare_matrix(Integer_Matrix<kind>& A,unsigned int n)
  {
    if (n >= A.nrow) throw std::invalid_argument("The row number argument is illegal for this matrix!");
    unsigned int i,j,v1,v2;
    bool jump;
    kind f;

    do {
      move_minimum(A,n);
      pivot_row(A,n,n);
      jump = false;
      for(i=n+1; i<A.nrow; ++i) {
        for(j=0; j<A.elements[i].size(); ++j) {
          if (A.elements[i][j].second == n) {
            jump = true;
            break;
          }
        }
        if (jump) break;
      }
      if (jump) continue;
      pivot_column(A,n,n);
      jump = false;
      for(i=0; i<A.elements[n].size(); ++i) {
        if (A.elements[n][i].second > n) {
          jump = true;
          break;
        }
      }
      if (jump) continue;
      if (A.divisible(n,&v1,&v2,&f)) {
        combine(A,v1,n,Integer_Matrix<kind>::unity,'r');
        combine(A,n,v2,-f,'c');
      }
      else {
        break;
      }
    } while(true);
  }

  template<class kind>
  unsigned int normalize(Integer_Matrix<kind>& A)
  {
    if (A.normalized) return 0;

    unsigned int i,j,p = 0,q = 0;
    bool quit;
    kind bt = Integer_Matrix<kind>::zero;

    do {
      p += 1;
      prepare_matrix(A,p-1);
      bt = A.get_diagonal(p-1);
      if (bt < Integer_Matrix<kind>::zero) {
        invert(A,p-1,'r');
        bt *= Integer_Matrix<kind>::neg1;
      }
      if (bt == Integer_Matrix<kind>::unity) q += 1;
      quit = true;
      for(i=p; i<A.nrow; ++i) {
        for(j=0; j<A.elements[i].size(); ++j) {
          if (A.elements[i][j].second >= p) {
            quit = false;
            break;
          }
        }
        if (!quit) break;
      }
      if (quit) break;
    } while(true);
    A.normalized = true;
    return q;
  }
}

