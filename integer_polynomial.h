#include "random.h"
#include "rational.h"

#ifndef _integerpolynomialh
#define _integerpolynomialh

namespace SYNARMOSMA {
  template<class kind>
  class Integer_Polynomial;

  template<class kind>
  std::ostream& operator <<(std::ostream&,const Integer_Polynomial<kind>&);

  template<class kind>
  bool operator ==(const Integer_Polynomial<kind>&,const Integer_Polynomial<kind>&);

  template<class kind>
  Integer_Polynomial<kind> operator +(const Integer_Polynomial<kind>&,const Integer_Polynomial<kind>&);

  template<class kind>
  Integer_Polynomial<kind> operator -(const Integer_Polynomial<kind>&,const Integer_Polynomial<kind>&);

  template<class kind>
  Integer_Polynomial<kind> operator *(kind,const Integer_Polynomial<kind>&);

  template<class kind>
  Integer_Polynomial<kind> operator *(const Integer_Polynomial<kind>&,const Integer_Polynomial<kind>&);

  template<class kind>
  // A class intended for low-degree polynomials since we store every coefficient, zero or not.
  class Integer_Polynomial {
   protected:
    unsigned int degree = 0;
    unsigned int characteristic = 0;
    bool irreducible = false;
    bool homogeneous = false;
    bool normed = false;
    std::vector<kind> terms;
  
    void initialize();
    inline void simplify();
    int write_terms(std::ofstream&) const;
    int read_terms(std::ifstream&);
   public:
    Integer_Polynomial();
    Integer_Polynomial(unsigned int);
    Integer_Polynomial(unsigned int,unsigned int);
    Integer_Polynomial(const std::vector<kind>&);
    Integer_Polynomial(const std::vector<kind>&,unsigned int);
    ~Integer_Polynomial();
    Integer_Polynomial& operator =(const Integer_Polynomial&);
    Integer_Polynomial& operator -(const Integer_Polynomial&);
    Integer_Polynomial(const Integer_Polynomial&);
    Integer_Polynomial<unsigned int> reduce(unsigned int);
    kind evaluate(kind);
    void generate(unsigned int);
    inline unsigned int get_degree() const {return degree;};
    inline bool is_null() const;
    kind get_value(unsigned int) const;
    void set_value(kind,unsigned int);
    void clear();
    int serialize(std::ofstream&) const;
    int deserialize(std::ifstream&);
    Integer_Polynomial<kind> derivative() const;
    friend std::ostream& operator << <>(std::ostream&,const Integer_Polynomial<kind>&);
    friend bool operator == <>(const Integer_Polynomial<kind>&,const Integer_Polynomial<kind>&);
    friend Integer_Polynomial<kind> operator +<>(const Integer_Polynomial<kind>&,const Integer_Polynomial<kind>&);
    friend Integer_Polynomial<kind> operator *<>(kind,const Integer_Polynomial<kind>&);
    friend Integer_Polynomial<kind> operator *<>(const Integer_Polynomial<kind>&,const Integer_Polynomial<kind>&);
  };

  template<class kind>
  void Integer_Polynomial<kind>::simplify()
  {
    unsigned int i,d = 0;
    if (characteristic > 0) {
      for(i=0; i<=degree; ++i) {
        terms[i] = terms[i] % characteristic;
      }
    }
    for(i=0; i<=degree; ++i) {
      if (terms[i] == kind(0)) continue;
      d = i;
    }
    if (d < degree) {
      std::vector<kind> nterms;

      for(i=0; i<=d; ++i) {
        nterms.push_back(terms[i]);
      }
      terms = nterms;
      degree = d;
    }
    if (terms[degree] == kind(1)) normed = true;
    if (terms[0] == kind(0)) homogeneous = true;
  }

  template<class kind> 
  bool Integer_Polynomial<kind>::is_null() const
  {
    if (terms.empty()) return true;
    unsigned int i;
    for(i=0; i<=degree; ++i) {
      if (terms[i] != kind(0)) return false;
    }
    return true;
  }

  template<class kind>
  std::ostream& operator <<(std::ostream& s,const Integer_Polynomial<kind>& source)
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
          s << "-" << -source.terms[1] << "*x ";
        }
        else {
          s << "-" << -source.terms[source.degree] << "*x^" << source.degree << " "; 	
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
            s << "- " << -source.terms[i] << "*x^" << i << " ";
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
            s << "- " << -source.terms[i] << "*x ";
          }	    	
        }  	
      }   	
    }
    if (source.terms[0] != 0) {
      if (source.terms[0] > 0) {
        s << "+ " << source.terms[0];
      }
      else {
        s << "- " << -source.terms[0]; 		
      }
    } 	
    return s;
  }

  template<class kind>
  bool operator ==(const Integer_Polynomial<kind>& p1,const Integer_Polynomial<kind>& p2)
  {
    if (p1.degree != p2.degree) return false;
    unsigned int i;
    for(i=0; i<p1.degree; ++i) {
      if (p1.terms[i] != p2.terms[i]) return false;
    }
    return true;
  }

  template<class kind>
  Integer_Polynomial<kind> operator -(const Integer_Polynomial<kind>& p1,const Integer_Polynomial<kind>& p2)
  {
    unsigned int i,mu = std::min(p1.degree,p2.degree);
    std::vector<kind> new_terms;

    for(i=0; i<=mu; ++i) {
      new_terms.push_back(p1.terms[i] - p2.terms[i]);
    }
    if (p1.degree > p2.degree) {
      for(i=1+mu; i<=p1.degree; ++i) {
        new_terms.push_back(p1.terms[i]);
      }
    }
    else if (p2.degree > p1.degree) {
      for(i=1+mu; i<=p2.degree; ++i) {
        new_terms.push_back(-p2.terms[i]);
      }
    }
    Integer_Polynomial<kind> output(new_terms);
    output.simplify();
    return output;
  }

  template<class kind>
  Integer_Polynomial<kind> operator +(const Integer_Polynomial<kind>& p1,const Integer_Polynomial<kind>& p2)
  {
    unsigned int i,mu = std::min(p1.degree,p2.degree);
    std::vector<kind> new_terms;

    for(i=0; i<=mu; ++i) {
      new_terms.push_back(p1.terms[i] + p2.terms[i]);
    }
    if (p1.degree > p2.degree) {
      for(i=1+mu; i<=p1.degree; ++i) {
        new_terms.push_back(p1.terms[i]);
      }
    }
    else if (p2.degree > p1.degree) {
      for(i=1+mu; i<=p2.degree; ++i) {
        new_terms.push_back(p2.terms[i]);
      }
    }
    Integer_Polynomial<kind> output(new_terms);
    output.simplify();
    return output;
  }

  template<class kind>
  Integer_Polynomial<kind> operator *(kind alpha,const Integer_Polynomial<kind>& p)
  {
    unsigned int i;
    std::vector<kind> new_terms;

    for(i=0; i<=p.degree; ++i) {
      new_terms.push_back(alpha*p.terms[i]);
    }
    Integer_Polynomial<kind> output(new_terms);
    output.simplify();
    return output;
  }

  template<class kind>
  Integer_Polynomial<kind> operator *(const Integer_Polynomial<kind>& p1,const Integer_Polynomial<kind>& p2)  
  {
    unsigned int i,j,k,mdegree = p1.degree + p2.degree;
    kind sum;
    std::vector<kind> new_terms;

    new_terms.push_back(p1.terms[0]*p2.terms[0]);
    for(i=0; i<mdegree-1; ++i) {
      // Need all the two factor additive partitions of "i"
      sum = 0;
      for(j=0; j<=p1.degree; ++j) {
        for(k=0; k<=p2.degree; ++k) {
          if ((j+k) == i) sum = sum + p1.terms[j]*p2.terms[k];
        }
      }
      new_terms.push_back(sum);
    }    
    new_terms.push_back(p1.terms[p1.degree]*p2.terms[p2.degree]);        
    Integer_Polynomial<kind> output(new_terms);
    output.simplify();
    return output;
  } 
} 
#endif





