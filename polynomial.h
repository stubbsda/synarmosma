#include "rational.h"

#ifndef _polynomialh
#define _polynomialh

template<class kind>
class Polynomial;

template<class kind>
std::ostream& operator <<(std::ostream&,const Polynomial<kind>&);

template<class kind>
Polynomial<kind> operator +(const Polynomial<kind>&,const Polynomial<kind>&);

template<class kind>
Polynomial<kind> operator -(const Polynomial<kind>&,const Polynomial<kind>&);

template<class kind>
Polynomial<kind> operator *(const Polynomial<kind>&,const Polynomial<kind>&);

template<class kind>
class Polynomial {
 private:
  unsigned int degree;
  std::vector<kind> terms;
  bool irreducible;
  bool homogeneous;
  bool normed;
  unsigned int characteristic;
  
  void initialize();
 public:
  Polynomial();
  Polynomial(unsigned int);
  Polynomial(unsigned int,unsigned int);
  Polynomial(const std::vector<kind>&);
  ~Polynomial();
  Polynomial& operator =(const Polynomial&);
  Polynomial& operator -(const Polynomial&);
  Polynomial(const Polynomial&);
  Polynomial<unsigned int> reduce(unsigned int);
  kind evaluate(kind);
  void factorize(std::vector<Polynomial>&);
  kind get_value(unsigned int) const;
  void set_value(kind,unsigned int);
  Polynomial<kind> derivative() const;
  friend std::ostream& operator << <>(std::ostream&,const Polynomial<kind>&);
  friend Polynomial<kind> operator +<>(const Polynomial<kind>&,const Polynomial<kind>&);
  friend Polynomial<kind> operator *<>(const Polynomial<kind>&,const Polynomial<kind>&);
};

template<class kind>
std::ostream& operator <<(std::ostream& s,const Polynomial<kind>& source)
{
  unsigned int i;

  if (source.terms[source.degree] > 0) {
    if (source.terms[source.degree] == 1) {
      if (source.degree == 1) {
        s << "x ";
      }
      else {
        s << "x^" << source.degree << " ";
      }
    }
    else {
      if (source.degree == 1) {
        s << source.terms[1] << "*x ";
      }
      else {
        s << source.terms[source.degree] << "*x^" << source.degree << " ";
      }
    }
  }
  else {
    if (source.terms[source.degree] == -1) {
      if (source.degree == 1) {
        s << "-x ";
      }
      else {
        s << "-x^" << source.degree << " ";
      }
    }
    else {
      if (source.degree == 1) {
        s << "-" << abs(source.terms[1]) << "*x ";
      }
      else {
        s << "-" << abs(source.terms[source.degree]) << "*x^" << source.degree << " "; 	
      }
    }
  }
  for(i=source.degree-1; i>=1; --i) {
    if (source.terms[i] == 0) continue;
    if (i > 1) {
      if (source.terms[i] == 1) {
        s << "x^" << i << " ";
      }
      else if (source.terms[i] == -1) {
        s << "-x^" << i << " ";
      }
      else {
        if (source.terms[i] > 0) {
          s << "+ " << source.terms[i] << "*x^" << i << " ";
        }
        else {
          s << "- " << abs(source.terms[i]) << "*x^" << i << " ";
        }	    	
      }
    }	
    else if (i == 1) {
      if (source.terms[i] == 1) {
        s << "+ x ";
      }
      else if (source.terms[i] == -1) {
        s << "- x ";
      }
      else {
        if (source.terms[i] > 0) {
          s << "+ " << source.terms[i] << "*x ";
        }
        else {
          s << "- " << abs(source.terms[i]) << "*x ";
        }	    	
      }  	
    }  	
  }
  if (source.terms[0] != 0) {
    if (source.terms[0] > 0) {
      s << "+ " << source.terms[0];
    }
    else {
      s << "- " << abs(source.terms[0]); 		
    }
  } 	
  return s;
}  
#endif





