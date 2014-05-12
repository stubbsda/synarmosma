#include "matrix.h"

template<>
const int Matrix<int>::zero = 0;
template<>
const int Matrix<int>::neg1 = -1;
template<>
const int Matrix<int>::unity = 1; 

template<>
const NTL::ZZ Matrix<NTL::ZZ>::zero = NTL::to_ZZ(0);
template<>
const NTL::ZZ Matrix<NTL::ZZ>::neg1 = NTL::to_ZZ(-1);
template<>
const NTL::ZZ Matrix<NTL::ZZ>::unity = NTL::to_ZZ(1); 

template<>
const double Matrix<double>::zero = 0.0;
template<>
const double Matrix<double>::neg1 = -1.0;
template<>
const double Matrix<double>::unity = 1.0; 

template<class kind>
Matrix<kind>::Matrix()
{
  normalized = false;

  nrow = 0;
  ncolumn = 0;
}

template<class kind>
Matrix<kind>::Matrix(unsigned int n)
{
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
Matrix<kind>::Matrix(unsigned int n,unsigned int m)
{
  normalized = false;
  nrow = n;
  ncolumn = m;
  elements = new std::vector<std::pair<kind,unsigned int> >[nrow];
}

template<class kind>
Matrix<kind>::Matrix(const Matrix<kind>& source)
{
  unsigned int i;
  delete[] elements;
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
  delete[] elements;
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
  delete[] elements;
}

template<class kind>
void Matrix<kind>::display() const
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
        std::cout << "%" << w << " ";
      }
      else {
        std::cout << w << " ";
      }
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}

template<>
bool Matrix<double>::divisible(unsigned int n,unsigned int* out1,unsigned int* out2,double* f) const
{
  // A meaningless operation for this kind of matrix...
  return true;
}


template<class kind>
bool Matrix<kind>::divisible(unsigned int n,unsigned int* out1,unsigned int* out2,kind* f) const
{
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
void Matrix<kind>::clear(bool deep)
{
  normalized = false;
  if (deep) {
    nrow = 0;
    ncolumn = 0;
    delete[] elements;
  }
  else {
    unsigned int i;
    for(i=0; i<nrow; ++i) {
      elements[i].clear();
    }
  }
}

template<class kind>
void Matrix<kind>::initialize(unsigned int n,unsigned int m)
{
  if (nrow == n) {
    clear(false);
  }
  else {
    clear(true);
  }
  nrow = n;
  ncolumn = m;
  elements = new std::vector<std::pair<kind,unsigned int> >[nrow];
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
  assert(A.ncolumn == b.size());

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
Matrix<kind>& operator *(const Matrix<kind>& A,const Matrix<kind>& B)
{
  assert(A.ncolumn == B.nrow);

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
void Matrix<kind>::set(unsigned int n,unsigned int m,kind v)
{
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
    elements[n][q].first = v;
  }
  else {
    elements[n].push_back(std::pair<kind,unsigned int>(v,m));
  }
}

template<class kind>
kind Matrix<kind>::get(unsigned int n,unsigned int m) const
{
  unsigned int i;
  kind output = zero;
  if (elements[n].empty()) return output;
  for(i=0; i<elements[n].size(); ++i) {
    if (elements[n][i].second == m) {
      output = elements[n][i].first;
      break;
    } 
  }
  return output;
}

template<class kind>
bool Matrix<kind>::empty_row(unsigned int n) const
{
  return elements[n].empty();
}

template<class kind>
kind Matrix<kind>::get_first_nonzero(unsigned int n) const
{
  if (elements[n].empty()) return zero;
  return elements[n][0].first;
}

template<class kind>
void permute(Matrix<kind>& A,unsigned int n,unsigned int m,char type)
{
  // If (type == c), then column(n) <-> column(m)
  // If (type == r), then row(n) <-> row(m)
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
void invert(Matrix<kind>& A,unsigned int n,char type)
{
  // If (type == c), then column(n) -> -1*column(n)
  // If (type == r), then row(n) -> -1*row(n)
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
void combine(Matrix<kind>& A,unsigned int n,unsigned int m,const kind& q,char type)
{
  // If (type == r), then row(n) -> row(n) + q*row(m)
  // If (type == c), then column(n) -> column(n) + q*column(m)
  unsigned int i,j,l = 0;
  bool found,fd1,fd2;
  kind wv;

  wv = 0;

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
void pivot_row(Matrix<kind>& A,unsigned int k,unsigned int l)
{
  unsigned int i,j;
  kind n,d,q;
  bool found;

  d = 0;

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
void pivot_column(Matrix<kind>& A,unsigned int k,unsigned int l)
{
  unsigned int i,j;
  kind n,d,q;
  bool found;

  d = 0;
 
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
void move_minimum(Matrix<kind>& A,unsigned int n)
{
  unsigned int i,j,im,jm,istart = 0,ulimit = A.nrow - n;  
  kind alpha,beta;
  bool full[ulimit];
  std::pair<kind,unsigned int>* VY = new std::pair<kind,unsigned int>[ulimit];
  std::vector<std::pair<kind,unsigned int> > VX;

  alpha = 0;

  for(i=n; i<A.nrow; ++i) {
    full[i-n] = false;
    for(j=0; j<A.elements[i].size(); ++j) {
      if (A.elements[i][j].second >= n) {
        VX.push_back(A.elements[i][j]);
      }
    }
    if (VX.empty()) continue;
    full[i-n] = true;
    alpha = abs(VX[0].first); jm = VX[0].second;
    for(j=1; j<VX.size(); ++j) {
      beta = abs(VX[j].first);
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
void prepare_matrix(Matrix<kind>& A,unsigned int n)
{
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
  kind bt;
  bool quit;

  bt = 0;

  do {
    p += 1;
    prepare_matrix(A,p-1);
    for(i=0; i<A.elements[p-1].size(); ++i) {
      if (A.elements[p-1][i].second == p-1) {
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
