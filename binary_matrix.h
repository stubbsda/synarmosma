#include "global.h"

#ifndef _bmatrixh
#define _bmatrixh

namespace SYNARMOSMA {
  class Binary_Matrix {
   private:
    std::vector<unsigned int>* elements;
    unsigned int nrow,ncolumn;

     void clear();
   public:
    Binary_Matrix();
    Binary_Matrix(unsigned int);
    Binary_Matrix(unsigned int,unsigned int);
    Binary_Matrix(const Binary_Matrix&);
    Binary_Matrix& operator =(const Binary_Matrix&);
    ~Binary_Matrix();
    int serialize(std::ofstream&) const;
    int deserialize(std::ifstream&);
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

