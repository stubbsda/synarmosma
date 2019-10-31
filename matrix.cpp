#include "matrix.h"

using namespace SYNARMOSMA;

template<class kind>
Matrix<kind>::Matrix()
{

}

template<class kind>
Matrix<kind>::Matrix(unsigned int n,bool identity)
{
  if (n == 0) throw std::invalid_argument("The matrix dimension must be greater than zero!");
  nrow = n;
  ncolumn = n;
  elements = new std::vector<std::pair<kind,unsigned int> >[nrow];
  if (identity) { 
    unsigned int i;
    for(i=0; i<nrow; ++i) {
      elements[i].push_back(std::pair<kind,unsigned int>(kind(1.0),i));
    }
  }
}

template<class kind>
Matrix<kind>::Matrix(unsigned int n,unsigned int m)
{
  if (n == 0 || m == 0) throw std::invalid_argument("The matrix dimensions must be greater than zero!");
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
      w = kind(0.0);
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
  if (mtype != 'r' && mtype != 'c') throw std::invalid_argument("The conversion type for the matrix must be row-oriented or column-oriented!");
  // This method supposes that A has already been allocated the appropriate amount
  // of memory...
  unsigned int i,j;

  for(i=0; i<nrow; ++i) {
    for(j=0; j<ncolumn; ++j) {
      A[nrow*i+j] = kind(0.0);
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
void Matrix<kind>::homotopy_scaling(kind a1,kind a2,Matrix<kind>* output) const
{
  // A method to compute output = a1*A + a2*I...
  if (output->nrow != nrow || output->ncolumn != ncolumn) std::invalid_argument("The output matrix matrix dimensions are illegal!");
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
void Matrix<kind>::multiply(const std::vector<kind>& b,std::vector<kind>& output) const
{
  if (ncolumn != b.size()) throw std::invalid_argument("The vector's length must equal the number of matrix columns!");
  unsigned int i,j;
  kind sum;

  output.clear();

  for(i=0; i<nrow; ++i) {
    sum = kind(0.0);
    for(j=0; j<elements[i].size(); ++j) {
      sum += elements[i][j].first*b[elements[i][j].second];
    }
    output.push_back(sum);
  }
}

template<class kind>
kind Matrix<kind>::determinant() const
{
  if (nrow != ncolumn) throw std::runtime_error("The determinant can only be computed for square matrices!");
  unsigned int i,j,k,q;
  bool found;
  std::vector<int> symbols,rest;
  std::vector<std::vector<int> > plist;
  kind product,output = kind(0.0);

  for(i=0; i<nrow; ++i) {
    rest.push_back(i);
  }
  permute(plist,symbols,rest,nrow);

  for(i=0; i<plist.size(); ++i) {
    product = (parity(plist[i]) == 1) ? kind(1.0) : kind(-1.0);
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
        product = kind(0.0);
        break;
      }
    }
    output += product;
  }

  return output;
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
bool Matrix<kind>::diagonally_dominant(unsigned int r) const 
{
  if (r >= nrow) throw std::invalid_argument("The row number argument is illegal for this matrix!");
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
bool Matrix<kind>::diagonally_dominant(unsigned int r,unsigned int c) const 
{
  if (r >= nrow) throw std::invalid_argument("The row number argument is illegal for this matrix!");
  if (c >= ncolumn) throw std::invalid_argument("The column number argument is illegal for this matrix!");
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
void Matrix<kind>::get_diagonal(std::vector<kind>& output) const
{
  unsigned int i,j;
  kind d;

  output.clear();

  for(i=0; i<nrow; ++i) {
    d = kind(0.0);
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
void Matrix<kind>::gauss_seidel_solver(std::vector<kind>& soln,const std::vector<kind>& source,double threshold,int max_its)
{
  if (soln.size() != nrow || source.size() != nrow) throw std::invalid_argument("The vector dimensions are wrong for the Gauss-Seidel solver!");
  unsigned int i,j;
  int its = 0;
  kind sum;
  double rho,error,error_new;
  std::vector<unsigned int> permutation;
  std::vector<kind> x_new,diagonal,x = soln,b = source;

  for(i=0; i<nrow; ++i) {
    x_new.push_back(kind(0.0));
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
    if (std::abs(diagonal[i]) < std::numeric_limits<double>::epsilon()) throw std::runtime_error("Diagonal element is zero!");
  }

  error = 0.0;
  for(i=0; i<nrow; ++i) {
    sum = kind(0.0);
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
      sum = kind(0.0);
      for(j=0; j<elements[i].size(); ++j) {
        if (elements[i][j].second == i) continue;
        sum += elements[i][j].first*x[elements[i][j].second];
      }
      x_new[i] = (b[i] - sum)/diagonal[i];
    }
    error_new = 0.0;
    for(i=0; i<nrow; ++i) {
      sum = kind(0.0);
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
    if (error_new > 10000.0) throw std::runtime_error("Iterations are diverging!");
    if (error_new > 2.0*error && error_new > 5.0) throw std::runtime_error("Iterations are diverging!");
    error = error_new;
    x = x_new;
    if (its >= max_its) break;
  } while(true);

  if (its == max_its && !(error < threshold)) throw std::runtime_error("Failed to converge after maximum iterations!");

  soln = x_new;
}

template<class kind>
void Matrix<kind>::clear(bool deep)
{
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
bool Matrix<kind>::consistent() const
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
unsigned int Matrix<kind>::eliminate_zeros()
{
  // This method finds ane eliminates zeros from the matrix elements...
  unsigned int i,j,nmodify = 0;
  std::vector<std::pair<kind,unsigned int> > rvector;

  for(i=0; i<nrow; ++i) {
    rvector.clear();
    for(j=0; j<elements[i].size(); ++j) {
      if (std::abs(elements[i][j].first) < std::numeric_limits<double>::epsilon()) {
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
kind Matrix<kind>::get_first_nonzero(unsigned int r) const
{
  if (r >= nrow) throw std::invalid_argument("The row number argument is illegal for this matrix!");
  if (elements[r].empty()) return 0.0;
  unsigned int i,min_index = ncolumn;
  kind output = 0.0;

  for(i=0; i<elements[r].size(); ++i) {
    if (elements[r][i].second < min_index && std::abs(elements[r][i].first) > std::numeric_limits<double>::epsilon()) {
      output = elements[r][i].first;
      min_index = elements[r][i].second;
    }
  }
  return output;
}

template<class kind>
unsigned int Matrix<kind>::number_nonzero() const
{
  unsigned int i,j,nzero = 0;

  for(i=0; i<nrow; ++i) {
    for(j=0; j<elements[i].size(); ++j) {
      if (std::abs(elements[i][j].first) > std::numeric_limits<double>::epsilon()) nzero++;
    }
  }

  return nzero;
}

template<class kind>
bool Matrix<kind>::empty_row(unsigned int r) const 
{
  if (r >= nrow) throw std::invalid_argument("The row number argument is illegal for this matrix!");
  unsigned int i;
  bool output = true;

  for(i=0; i<elements[r].size(); ++i) {
    if (std::abs(elements[r][i].first) > std::numeric_limits<double>::epsilon()) {
      output = false;
      break;
    }
  }
  return output;
}

template<class kind>
void Matrix<kind>::initialize(unsigned int n,unsigned int m)
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
  if (A.ncolumn != b.size()) throw std::invalid_argument("The matrix and vector dimensions don't conform for multiplication!");
  unsigned int i,j;
  kind sum;
  std::vector<kind> output;

  for(i=0; i<A.nrow; ++i) {
    sum = kind(0.0);
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
  if (A.nrow != B.nrow || A.ncolumn != B.ncolumn) throw std::invalid_argument("The dimensions of the two matrices don't conform for addition!");
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
  output.eliminate_zeros();
  return output;
}

template<class kind>
Matrix<kind>& operator *(const Matrix<kind>& A,const Matrix<kind>& B)
{
  if (A.ncolumn != B.nrow) throw std::invalid_argument("The dimensions of the two matrices don't conform for multiplication!");
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
      sum = kind(0.0);
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
void Matrix<kind>::get_row(std::vector<kind>& v,std::vector<unsigned int>& c,unsigned int r) const
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
unsigned int Matrix<kind>::optimize_dominance(std::vector<unsigned int>& shift)
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
  unsigned int ndominant = 0;
  for(i=0; i<nrow; ++i) {
    if (diagonally_dominant(i)) ndominant++;
  }
  return ndominant;
}
