#include "matrix.h"

using namespace SYNARMOSMA;

extern Random RND;

template<class kind>
Matrix<kind>::Matrix()
{
  normalized = false;

  nrow = 0;
  ncolumn = 0;
}

template<class kind>
Matrix<kind>::Matrix(int n)
{
  assert(n > 0);
  unsigned int i;

  normalized = false;
  nrow = n;
  ncolumn = n;
  elements = new std::vector<std::pair<kind,unsigned int> >[nrow];
  for(i=0; i<nrow; ++i) {
    elements[i].push_back(std::pair<kind,unsigned int>(unity,i));
  }
}

template<class kind>
Matrix<kind>::Matrix(int n,int m)
{
  assert(n > 0 && m > 0);
  normalized = false;
  nrow = n;
  ncolumn = m;
  elements = new std::vector<std::pair<kind,unsigned int> >[nrow];
}

template<class kind>
Matrix<kind>::Matrix(const Matrix<kind>& source)
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
Matrix<kind>& Matrix<kind>::operator =(const Matrix<kind>& source)
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
Matrix<kind>::~Matrix()
{
  if (nrow > 0) delete[] elements;
}

template<class kind>
void Matrix<kind>::display(std::ostream& s) const
{
  unsigned int i,j,k;
  kind w;
  bool found;

  for(i=0; i<nrow; ++i) {
    for(j=0; j<ncolumn; ++j) {
      w = zero;
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
void Matrix<kind>::convert(kind* A,char mtype) const
{
#ifdef DEBUG
  assert(mtype == 'r' || mtype == 'c');
#endif
  // This method supposes that A has already been allocated the appropriate amount
  // of memory...
  unsigned int i,j;

  for(i=0; i<nrow; ++i) {
    for(j=0; j<ncolumn; ++j) {
      A[nrow*i+j] = zero;
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

namespace SYNARMOSMA {
  // A divisibility test is meaningless in the context of fields like \mathbb{R} and \mathbb{C}...
  template<>
  bool Matrix<double>::divisible(int n,unsigned int* out1,unsigned int* out2,double* f) const
  {
    assert(n >= 0);
    return true;
  }

  template<>
  bool Matrix<std::complex<double> >::divisible(int n,unsigned int* out1,unsigned int* out2,std::complex<double>* f) const
  {
    assert(n >= 0);
    return true;
  }
}

template<class kind>
bool Matrix<kind>::divisible(int n,unsigned int* out1,unsigned int* out2,kind* f) const
{
  assert(n >= 0);
  unsigned int i,j;
  kind b1,b2 = Matrix<kind>::zero;

  for(i=0; i<elements[n].size(); ++i) {
    if (elements[n][i].second == n) b2 = elements[n][i].first;
  }
  if (b2 == Matrix<kind>::zero) return false;
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
void Matrix<kind>::homotopy_scaling(kind a1,kind a2,Matrix<kind>* output) const
{
  // A method to compute output = a1*A + a2*I...
#ifdef DEBUG
  assert(output->nrow == nrow && output->ncolumn == ncolumn);
#endif
  unsigned int i,j;
  bool diagonal;

  output->clear(false);
  for(i=0; i<nrow; ++i) {
    diagonal = false;
    for(j=0; j<elements[i].size(); ++j) {
      if (elements[i][j].second == i) {
        diagonal = true;
        output->set(i,i,a1*elements[i][j].first + a2);
        continue;
      }
      output->set(i,elements[i][j].second,a1*elements[i][j].first);
    }
    if (!diagonal) output->set(i,i,a2);
  }
}

template<class kind>
void Matrix<kind>::increment(const Matrix<kind>& arg)
{
#ifdef DEBUG
  assert(nrow == arg.nrow && ncolumn == arg.ncolumn);
#endif
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
void Matrix<kind>::multiply(const std::vector<kind>& b,std::vector<kind>& output) const
{
#ifdef DEBUG
  assert(ncolumn == b.size());
#endif
  unsigned int i,j;
  kind sum;

  output.clear();

  for(i=0; i<nrow; ++i) {
    sum = Matrix<kind>::zero;
    for(j=0; j<elements[i].size(); ++j) {
      sum += elements[i][j].first*b[elements[i][j].second];
    }
    output.push_back(sum);
  }
}

template<class kind>
int Matrix<kind>::number_nonzero() const
{
  unsigned int i,nzero = 0;

  for(i=0; i<nrow; ++i) {
    nzero += elements[i].size();
  }
  return nzero;
}

template<class kind>
kind Matrix<kind>::determinant() const
{
  assert(nrow == ncolumn);
  unsigned int i,j,k,q;
  bool found;
  std::vector<int> symbols,rest;
  std::vector<std::vector<int> > plist;
  kind product,output = zero;

  for(i=0; i<nrow; ++i) {
    rest.push_back(i);
  }
  permute(plist,symbols,rest,nrow);

  for(i=0; i<plist.size(); ++i) {
    product = (parity(plist[i]) == 1) ? unity : neg1;
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
        product = zero;
        break;
      }
    }
    output += product;
  }

  return output;
}

namespace SYNARMOSMA {
  template<>
  bool Matrix<NTL::ZZ>::diagonally_dominant() const
  {
    unsigned int i,j;
    double sum,d,q;
    bool output = true;

#ifndef VERBOSE
    unsigned int nf = 0;
    for(i=0; i<nrow; ++i) {
      d = 0.0;
      sum = 0.0;
      for(j=0; j<elements[i].size(); ++j) {
        NTL::conv(elements[i][j].first,q);
        q = std::abs(q);
        if (elements[i][j].second == i) {
          d = q;
          continue;
        }
        sum += q;
      }
      if (d < sum || std::abs(d) < std::numeric_limits<double>::epsilon()) {
        nf++;
        output = false;
      }
      else {
        std::cout << "DOMINANT ";
      }
      std::cout << 1+i << "  " << d << "  " << sum << "   " << sum/d << std::endl;
    }
    std::cout << "There are " << nrow - nf << " diagonally dominant rows." << std::endl;
#else
    for(i=0; i<nrow; ++i) {
      d = 0.0;
      sum = 0.0;
      for(j=0; j<elements[i].size(); ++j) {
        NTL::conv(elements[i][j].first,q);
        q = std::abs(q);
        if (elements[i][j].second == i) {
          d = q;
          continue;
        }
        sum += q;
      }
      if (d < sum || std::abs(d) < std::numeric_limits<double>::epsilon()) return false;
    }
#endif
    return output;
  }

  template<>
  bool Matrix<NTL::ZZ>::diagonally_dominant(int r) const 
  {
    assert(r >= 0);
    unsigned int i;
    double q,d = 0.0,sum = 0.0;

    for(i=0; i<elements[r].size(); ++i) {
      NTL::conv(elements[r][i].first,q);
      q = std::abs(q);
      if (elements[r][i].second == r) {
        d = q;
        continue;
      }
      sum += q;
    }
    if (d < sum || std::abs(d) < std::numeric_limits<double>::epsilon()) return false;

    return true; 
  }

  template<>
  bool Matrix<NTL::ZZ>::diagonally_dominant(int r,int c) const 
  {
    assert(r >= 0 && c >= 0);
    unsigned int i;
    double q,d = 0.0,sum = 0.0;

    for(i=0; i<elements[r].size(); ++i) {
      NTL::conv(elements[r][i].first,q);
      q = std::abs(q);
      if (elements[r][i].second == c) {
        d = q;
        continue;
      }
      sum += q;
    }
    if (d < sum || std::abs(d) < std::numeric_limits<double>::epsilon()) return false;

    return true; 
  }
}

template<class kind>
bool Matrix<kind>::diagonally_dominant() const
{
  unsigned int i,j;
  double sum,d,q;
  bool output = true;

#ifndef VERBOSE
  unsigned int nf = 0;
  for(i=0; i<nrow; ++i) {
    d = 0.0;
    sum = 0.0;
    for(j=0; j<elements[i].size(); ++j) {
      q = std::abs(elements[i][j].first);
      if (elements[i][j].second == i) {
        d = q;
        continue;
      }
      sum += q;
    }
    if (d < sum || std::abs(d) < std::numeric_limits<double>::epsilon()) {
      nf++; 
      output = false;
    }
    else {
      std::cout << "DOMINANT ";
    }
    std::cout << 1+i << "  " << d << "  " << sum << "   " << sum/d << std::endl;
  }
  std::cout << "There are " << nrow - nf << " diagonally dominant rows." << std::endl;
#else
  for(i=0; i<nrow; ++i) {
    d = 0.0;
    sum = 0.0;
    for(j=0; j<elements[i].size(); ++j) {
      q = std::abs(elements[i][j].first);
      if (elements[i][j].second == i) {
        d = q;
        continue;
      }
      sum += q;
    }
    if (d < sum || std::abs(d) < std::numeric_limits<double>::epsilon()) return false;
  }
#endif
  return output;
}

template<class kind>
bool Matrix<kind>::diagonally_dominant(int r) const 
{
  assert(r >= 0);
  unsigned int i;
  double q,d = 0.0,sum = 0.0;

  for(i=0; i<elements[r].size(); ++i) {
    q = std::abs(elements[r][i].first);
    if (elements[r][i].second == r) {
      d = q;
      continue;
    }
    sum += q;
  }
  if (d < sum || std::abs(d) < std::numeric_limits<double>::epsilon()) return false;

  return true; 
}

template<class kind>
bool Matrix<kind>::diagonally_dominant(int r,int c) const 
{
  assert(r >= 0 && c >= 0);
  unsigned int i;
  double q,d = 0.0,sum = 0.0;

  for(i=0; i<elements[r].size(); ++i) {
    q = std::abs(elements[r][i].first);
    if (elements[r][i].second == c) {
      d = q;
      continue;
    }
    sum += q;
  }
  if (d < sum || std::abs(d) < std::numeric_limits<double>::epsilon()) return false;

  return true; 
}

template<class kind>
void Matrix<kind>::get_diagonal(std::vector<kind>& diag) const
{
  unsigned int i,j;
  kind d;

  diag.clear();

  for(i=0; i<nrow; ++i) {
    d = zero;
    for(j=0; j<elements[i].size(); ++j) {
      if (elements[i][j].second == i) {
        d = elements[i][j].first;
        break;
      }
    }
    diag.push_back(d);
  }
}

namespace SYNARMOSMA {
  // This method only makes sense in the context of a field, so we will create null versions for the int 
  // and NTL::ZZ instantiations of the class
  template<>
  int Matrix<int>::gauss_seidel_solver(std::vector<int>& x,const std::vector<int>& b,double threshold,int max_its)
  {
    return -1;
  }

  template<>
  int Matrix<NTL::ZZ>::gauss_seidel_solver(std::vector<NTL::ZZ>& x,const std::vector<NTL::ZZ>& b,double threshold,int max_its)
  {
    return -1;
  }
}

template<class kind>
int Matrix<kind>::gauss_seidel_solver(std::vector<kind>& soln,const std::vector<kind>& source,double threshold,int max_its)
{
#ifdef DEBUG
  assert(soln.size() == nrow);
  assert(source.size() == nrow);
#endif
  unsigned int i,j;
  int its = 0;
  kind sum;
  double rho,error,error_new;
  std::vector<unsigned int> permutation;
  std::vector<kind> x_new,diagonal,x = soln,b = source;

  for(i=0; i<nrow; ++i) {
    x_new.push_back(zero);
  }

  if (!diagonally_dominant(true)) {
    j = optimize_dominance(permutation);
    // Fix the ordering of the source vector and initial guess 
    // vectors...
    for(i=0; i<nrow; ++i) {
      b[permutation[i]] = source[i];
      x[permutation[i]] = soln[i];
    }
#ifdef VERBOSE
    if (j < nrow) std::cout << "Matrix is not diagonally dominant (" << j << "/" << nrow << "), convergence not guaranteed!" << std::endl;
#endif
  }

  get_diagonal(diagonal);
  for(i=0; i<nrow; ++i) {
#ifdef VERBOSE
    //std::cout << "Diagonal element of row " << i << " is " << diagonal[i] << std::endl;
#endif
    if (std::abs(diagonal[i]) < std::numeric_limits<double>::epsilon()) {
#ifdef VERBOSE
      std::cout << "Illegal diagonal element value, exiting!" << std::endl;
#endif
      return -1;
    }
  }

  error = 0.0;
  for(i=0; i<nrow; ++i) {
    sum = zero;
    for(j=0; j<elements[i].size(); ++j) {
      sum += elements[i][j].first*x[elements[i][j].second];
    }
    rho = std::abs(sum - b[i]);
    error += rho*rho;
  }
  error = std::sqrt(error);

#ifdef VERBOSE
  std::cout << "For iteration " << its << " the solution error of the linear system is " << error << std::endl;  
#endif

  do {
    its++;
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i,j,sum)
#endif
    for(i=0; i<nrow; ++i) {
      sum = zero;
      for(j=0; j<elements[i].size(); ++j) {
        if (elements[i][j].second == i) continue;
        sum += elements[i][j].first*x[elements[i][j].second];
      }
      x_new[i] = (b[i] - sum)/diagonal[i];
    }
    error_new = 0.0;
    for(i=0; i<nrow; ++i) {
      sum = zero;
      for(j=0; j<elements[i].size(); ++j) {
        sum += elements[i][j].first*x_new[elements[i][j].second];
      }
      rho = std::abs(sum - b[i]);
      error_new += rho*rho;
    }
    error_new = std::sqrt(error_new);
#ifdef VERBOSE
    std::cout << "For iteration " << its << " the solution error of the linear system is " << error_new << std::endl;
#endif
    if (error_new < threshold) break;
    if (error_new > 10000.0) return -1;
    if (error_new > 2.0*error && error_new > 5.0) return -1;
    error = error_new;
    x = x_new;
    if (its >= max_its) break;
  } while(true);

  if (its == max_its && !(error < threshold)) return -1;

  soln = x_new;

  return its;
}

template<class kind>
void Matrix<kind>::clear()
{
  if (nrow > 0) delete[] elements;
  nrow = 0;
  ncolumn = 0;
}

template<class kind>
void Matrix<kind>::clear(bool deep)
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

namespace SYNARMOSMA {
  template<>
  int Matrix<NTL::ZZ>::write_elements(std::ofstream& s) const
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
int Matrix<kind>::write_elements(std::ofstream& s) const
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
int Matrix<kind>::serialize(std::ofstream& s) const
{
  int count = 0;

  s.write((char*)(&nrow),sizeof(int)); count += sizeof(int);
  s.write((char*)(&ncolumn),sizeof(int)); count += sizeof(int);
  count += write_elements(s);

  return count;
}

namespace SYNARMOSMA {
  template<>
  int Matrix<NTL::ZZ>::read_elements(std::ifstream& s)
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
int Matrix<kind>::read_elements(std::ifstream& s)
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
int Matrix<kind>::deserialize(std::ifstream& s)
{
  int count = 0;

  clear();

  s.read((char*)(&nrow),sizeof(int)); count += sizeof(int);
  s.read((char*)(&ncolumn),sizeof(int)); count += sizeof(int);
  elements = new std::vector<std::pair<kind,unsigned int> >[nrow];
  count += read_elements(s);

  return count;
}

template<class kind>
void Matrix<kind>::initialize(int n,int m)
{
  assert(n > 0 && m > 0);
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
void Matrix<kind>::transpose(const Matrix<kind>& source)
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
std::vector<kind>& operator *(const Matrix<kind>& A,const std::vector<kind>& b)
{
#ifdef DEBUG
  assert(A.ncolumn == b.size());
#endif
  unsigned int i,j;
  kind sum;
  std::vector<kind> output;

  for(i=0; i<A.nrow; ++i) {
    sum = Matrix<kind>::zero;
    for(j=0; j<A.elements[i].size(); ++j) {
      sum += A.elements[i][j].first*b[A.elements[i][j].second];
    }
    output.push_back(sum);
  }
  return output;
}

template<class kind>
Matrix<kind>& operator +(const Matrix<kind>& A,const Matrix<kind>& B)
{
#ifdef DEBUG
  assert(A.nrow == B.nrow && A.ncolumn == B.ncolumn);
#endif
  unsigned int i,j,k,in1,cvalue;
  bool found;
  Matrix<kind> output(A.nrow,A.ncolumn);

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
  return output;
}

template<class kind>
Matrix<kind>& operator *(const Matrix<kind>& A,const Matrix<kind>& B)
{
#ifdef DEBUG
  assert(A.ncolumn == B.nrow);
#endif
  unsigned int i,j,k,l,in1;
  kind sum;
  Matrix<kind> output(A.nrow,B.ncolumn);

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
      sum = Matrix<kind>::zero;
      for(k=0; k<A.elements[i].size(); ++k) {
        in1 = A.elements[i][k].second;
        for(l=0; l<BC_row[j].size(); ++l) {
          if (BC_row[j][l] == in1) {
            sum += A.elements[i][k].first*BC[j][l];
          }
        }
      }
      if (sum != Matrix<kind>::zero) output.elements[i].push_back(std::pair<kind,unsigned int>(sum,j));
    }
  }
  delete[] BC;
  delete[] BC_row;
  return output;
}

template<class kind>
void Matrix<kind>::increment(int n,int m,kind v)
{
  assert(n >= 0 && m >= 0);
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
    elements[n][q].first += v;
  }
  else {
    elements[n].push_back(std::pair<kind,unsigned int>(v,m));
  }
}

template<class kind>
void Matrix<kind>::get_row(std::vector<kind>& v,std::vector<unsigned int>& c,int rnumber) const
{
  unsigned int i;
  v.clear();
  c.clear();
  for(i=0; i<elements[rnumber].size(); ++i) {
    v.push_back(elements[rnumber][i].first);
    c.push_back(elements[rnumber][i].second);
  }
}

template<class kind>
bool Matrix<kind>::empty_row(int n) const
{
  assert(n >= 0);
  return elements[n].empty();
}

template<class kind>
kind Matrix<kind>::get_first_nonzero(int n) const
{
  assert(n >= 0);
  if (elements[n].empty()) return zero;
  return elements[n][0].first;
}

template<class kind>
double Matrix<kind>::dispersion() const
{
  // Measures the percentage of rows which do not have a single dominant 
  // column
  unsigned int i,j,nc = 0;
  bool found;

  for(i=0; i<nrow; ++i) {
    found = false;
    for(j=0; j<ncolumn; ++j) {
      if (diagonally_dominant(i,j)) {
        found = true;
        break;
      }
    }
    if (!found) nc++;
  }
  return double(nc)/double(nrow);
}

template<class kind>
int Matrix<kind>::optimize_dominance(std::vector<unsigned int>& shift)
{
  // A method to permute the matrix rows so that the diagonal dominance of the 
  // matrix is improved
  unsigned int i,j,q,np,current,next,permutation[nrow],its = 0;
  bool found,output;
  std::vector<std::pair<kind,unsigned int> > temp;

  for(i=0; i<nrow; ++i) {
    permutation[i] = i;
  }

  do {
    np = 0;
    for(i=0; i<nrow; ++i) {
      if (diagonally_dominant(i)) continue;
      // This row is a problem, see if there is another row which can be swapped with 
      // it...
      found = false;
      for(j=0; j<i; ++j) {
        if (diagonally_dominant(j,i)) {
          temp = elements[i];
          elements[i] = elements[j];
          elements[j] = temp;
          q = permutation[i];
          permutation[i] = permutation[j];
          permutation[j] = q;
          found = true;
          np++;
        }
      }
      if (found) continue;
      for(j=1+i; j<nrow; ++j) {
        if (diagonally_dominant(j,i)) {
          temp = elements[i];
          elements[i] = elements[j];
          elements[j] = temp;
          q = permutation[i];
          permutation[i] = permutation[j];
          permutation[j] = q;
          found = true;
          np++;
        }      
      }
    }
    output = diagonally_dominant();
    its++;
    if (output || its > nrow || np == 0) break;
  } while(true);

  shift.clear();
  for(i=0; i<nrow; ++i) {
    current = i;
    do {
      next = permutation[current];
      if (next == i) {
        shift.push_back(current);
        break;
      }
      current = next;
    } while(true);
  }
  if (output) return nrow;
  // Compute the number of rows that are diagonally dominant and return 
  // this number
  int ndominant = 0;
  for(i=0; i<nrow; ++i) {
    if (diagonally_dominant(i)) ndominant++;
  }
  return ndominant;
}

namespace SYNARMOSMA {
  template<class kind>
  void permute(Matrix<kind>& A,int n,int m,char type)
  {
    // If (type == c), then column(n) <-> column(m)
    // If (type == r), then row(n) <-> row(m)
    assert(n >= 0 && m >= 0);
    unsigned int i,j;
    std::vector<std::pair<kind,unsigned int> > pr;

    if (type == 'c') {
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
      pr = A.elements[n];
      A.elements[n] = A.elements[m];
      A.elements[m] = pr;
    }
  }

  template<class kind>
  void invert(Matrix<kind>& A,int n,char type)
  {
    // If (type == c), then column(n) -> -1*column(n)
    // If (type == r), then row(n) -> -1*row(n)
    assert(n >= 0);
    unsigned int i,j;

    if (type == 'c') {
      for(i=0; i<A.nrow; ++i) {
        for(j=0; j<A.elements[i].size(); ++j) {
          if (A.elements[i][j].second == n) A.elements[i][j].first *= Matrix<kind>::neg1;
        }
      }
    }
    else {
      for(i=0; i<A.elements[n].size(); ++i) {
        A.elements[n][i].first *= Matrix<kind>::neg1;
      }
    }
  }

  template<class kind>
  void combine(Matrix<kind>& A,int n,int m,const kind& q,char type)
  {
    // If (type == r), then row(n) -> row(n) + q*row(m)
    // If (type == c), then column(n) -> column(n) + q*column(m)
    assert(n >= 0 && m >= 0);
    unsigned int i,j,l = 0;
    bool found,fd1,fd2;
    kind wv = Matrix<kind>::zero;

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
          if (A.elements[i][l].first == -wv) {
            A.elements[i].erase(A.elements[i].begin() + l);
          }
          else {
            A.elements[i][l].first += wv;
          }
        }
      }
    }
    else {
      for(i=0; i<A.elements[m].size(); ++i) {
        l = A.elements[m][i].second;
        found = false;
        for(j=0; j<A.elements[n].size(); ++j) {
          if (A.elements[n][j].second == l) {
            if (A.elements[n][j].first == -q*A.elements[m][i].first) {
              A.elements[n].erase(A.elements[n].begin() + j);
            }
            else {
              A.elements[n][j].first += q*A.elements[m][i].first;
            }
            found = true;
            break;
          }
        }
        if (!found) A.elements[n].push_back(std::pair<kind,unsigned int>(q*A.elements[m][i].first,l));
      }
    }
  }

  template<class kind>
  void pivot_row(Matrix<kind>& A,int k,int l)
  {
    assert(k >= 0 && l >= 0);
    unsigned int i,j;
    bool found;
    kind n,q,d = Matrix<kind>::zero;

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
  void pivot_column(Matrix<kind>& A,int k,int l)
  {
    assert(k >= 0 && l >= 0);
    unsigned int i,j;
    bool found;
    kind n,q,d = Matrix<kind>::zero;
 
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
  void move_minimum(Matrix<kind>& A,int n)
  {
    assert(n >= 0);
    unsigned int i,j,im,jm,istart = 0,ulimit = A.nrow - n;  
    kind beta,alpha = Matrix<kind>::zero;
    bool full[ulimit];
    std::pair<kind,unsigned int>* VY = new std::pair<kind,unsigned int>[ulimit];
    std::vector<std::pair<kind,unsigned int> > VX;

    for(i=n; i<A.nrow; ++i) {
      full[i-n] = false;
      for(j=0; j<A.elements[i].size(); ++j) {
        if (A.elements[i][j].second >= n) {
          VX.push_back(A.elements[i][j]);
        }
      }
      if (VX.empty()) continue;
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
  void prepare_matrix(Matrix<kind>& A,int n)
  {
    assert(n >= 0);
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
        combine(A,v1,n,Matrix<kind>::unity,'r');
        combine(A,n,v2,-f,'c');
      }
      else {
        break;
      }
    } while(true);
  }

  template<class kind>
  unsigned int normalize(Matrix<kind>& A)
  {
    if (A.normalized) return 0;

    unsigned int i,j,p = 0,q = 0;
    bool quit;
    kind bt = Matrix<kind>::zero;

    do {
      p += 1;
      prepare_matrix(A,p-1);
      for(i=0; i<A.elements[p-1].size(); ++i) {
        if (A.elements[p-1][i].second == (p - 1)) {
          bt = A.elements[p-1][i].first;
          break;
        }
      }
      if (bt < Matrix<kind>::zero) {
        invert(A,p-1,'r');
      }
      else if (bt == Matrix<kind>::unity) {
        q += 1;
      }
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

