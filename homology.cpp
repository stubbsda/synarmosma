#include "homology.h"

Homology::Homology()
{

}

Homology::Homology(FIELD f)
{
  field = f;
}

Homology::~Homology()
{

}

unsigned int Homology::normalize_operator(const std::vector<signed int>* boundary_operator,unsigned int r,unsigned int c,std::vector<unsigned int>& torsion) const
{
  unsigned int i,j,s,rank = 0;
  std::string cx;
  std::stringstream ss;

  if (field == INTEGER) {
    Matrix<int> A(r,c),Q(r,c),NQ(r,c),R(r,c),NR(r,c);
    int v;

    for(i=0; i<r; ++i) {
      for(j=0; j<boundary_operator[i].size(); ++j) {
        s = std::abs(boundary_operator[i][j]);
        v = (boundary_operator[i][j] > 0) ? 1 : -1;
        A.set(i,s-1,v);
      }
    }
    normalize(A,Q,NQ,R,NR);
    for(i=0; i<r; ++i) {
      if (A.empty_row(i)) break;
      s = (unsigned) A.get_first_nonzero(i);
      if (s == 1) {
        rank += 1;
      }
      else {
        torsion.push_back(s);
      }
    }
  }
  else if (field == ZZ) {
    Matrix<NTL::ZZ> A(r,c),Q(r,c),NQ(r,c),R(r,c),NR(r,c);
    NTL::ZZ v;

    for(i=0; i<r; ++i) {
      for(j=0; j<boundary_operator[i].size(); ++j) {
        s = std::abs(boundary_operator[i][j]);
        v = (boundary_operator[i][j] > 0) ? Matrix<NTL::ZZ>::unity : Matrix<NTL::ZZ>::neg1;
        A.set(i,s-1,v); 
      }
    }
    normalize(A,Q,NQ,R,NR);
    for(i=0; i<r; ++i) {
      if (A.empty_row(i)) break;
      v = A.get_first_nonzero(i);
      if (v == Matrix<NTL::ZZ>::unity) {
        rank += 1;
      }
      else {
        ss << v;
        cx = ss.str();
        torsion.push_back(boost::lexical_cast<unsigned int>(cx.c_str()));
        ss.str("");
      }
    }
  }
  else {
    // The Galois field GF2
    Binary_Matrix A(r,c);

    for(i=0; i<r; ++i) {
      for(j=0; j<boundary_operator[i].size(); ++j) {
        s = std::abs(boundary_operator[i][j]);
        if (s != 0) A.set(i,s-1);
      }
    }
    rank = A.rank();
  }
  return rank;
}

void Homology::betti_numbers(std::vector<unsigned int>& bnumbers) const
{
  unsigned int i;
  bnumbers.clear();
  for(i=0; i<sequence.size(); ++i) {
    bnumbers.push_back(sequence[i].get_rank());
  }
}

void Homology::append(const Group& g)
{
  sequence.push_back(g);
}

void Homology::clear()
{
  sequence.clear();
}


