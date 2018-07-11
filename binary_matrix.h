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
    Binary_Matrix(int);
    Binary_Matrix(int,int);
    Binary_Matrix(const Binary_Matrix&);
    Binary_Matrix& operator =(const Binary_Matrix&);
    ~Binary_Matrix();
    int serialize(std::ofstream&) const;
    int deserialize(std::ifstream&);
    void initialize(int,int);
    void get_row(int,bool*) const;
    bool get(int,int) const;
    void set(int,int);
    void unset(int,int);
    int rank() const;
    friend std::ostream& operator <<(std::ostream&,const Binary_Matrix&);
    friend Binary_Matrix operator *(const Binary_Matrix&,const Binary_Matrix&);
    friend Binary_Matrix operator +(const Binary_Matrix&,const Binary_Matrix&);
    friend class Graph;
  };
}
#endif

