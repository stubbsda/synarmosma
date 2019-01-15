#include "random.h"
#include "rational.h"
#include "integer_polynomial.h"

#ifndef _polynomialh
#define _polynomialh

namespace SYNARMOSMA {
  template<class kind>
  class Polynomial;

  template<class kind>
  std::ostream& operator <<(std::ostream&,const Polynomial<kind>&);

  template<class kind>
  bool operator ==(const Polynomial<kind>&,const Polynomial<kind>&);

  template<class kind>
  Polynomial<kind> operator +(const Polynomial<kind>&,const Polynomial<kind>&);

  template<class kind>
  Polynomial<kind> operator -(const Polynomial<kind>&,const Polynomial<kind>&);

  template<class kind>
  Polynomial<kind> operator *(kind,const Polynomial<kind>&);

  template<class kind>
  Polynomial<kind> operator *(const Polynomial<kind>&,const Polynomial<kind>&);

  template<class kind>
  // A class intended for low-degree polynomials since we store every coefficient, zero or not.
  class Polynomial {
   protected:
    unsigned int degree = 0;
    bool irreducible = false;
    bool homogeneous = false;
    bool monic = false;
    std::vector<kind> terms;
  
    void initialize();
    inline void simplify();
    int write_terms(std::ofstream&) const;
    int read_terms(std::ifstream&);
   public:
    Polynomial();
    Polynomial(unsigned int);
    Polynomial(unsigned int,unsigned int);
    Polynomial(const std::vector<kind>&);
    Polynomial(const std::vector<kind>&,unsigned int);
    ~Polynomial();
    Polynomial& operator =(const Polynomial&);
    Polynomial& operator -(const Polynomial&);
    Polynomial(const Polynomial&);
    Integer_Polynomial<int> reduce(unsigned int);
    kind evaluate(kind);
    void generate(unsigned int);
    inline unsigned int get_degree() const {return degree;};
    inline bool is_null() const;
    kind get_value(unsigned int) const;
    void set_value(kind,unsigned int);
    void clear();
    int serialize(std::ofstream&) const;
    int deserialize(std::ifstream&);
    Polynomial<kind> derivative() const;
    friend std::ostream& operator << <>(std::ostream&,const Polynomial<kind>&);
    friend bool operator == <>(const Polynomial<kind>&,const Polynomial<kind>&);
    friend Polynomial<kind> operator +<>(const Polynomial<kind>&,const Polynomial<kind>&);
    friend Polynomial<kind> operator *<>(kind,const Polynomial<kind>&);
    friend Polynomial<kind> operator *<>(const Polynomial<kind>&,const Polynomial<kind>&);
  };

  template<>
  void Polynomial<Rational>::simplify()
  {
    unsigned int i,d = 0;
    const Rational unity(1);

    for(i=0; i<=degree; ++i) {
      if (terms[i].is_null()) continue;
      d = i;
    }
    if (d < degree) {
      std::vector<Rational> nterms;

      for(i=0; i<=d; ++i) {
        nterms.push_back(terms[i]);
      }
      terms = nterms;
      degree = d;
    }
    if (terms[degree] == unity) monic = true;
    if (terms[0].is_null()) homogeneous = true;
  }

  template<class kind>
  void Polynomial<kind>::simplify()
  {
    unsigned int i,d = 0;
    const kind unity = 1.0;

    for(i=0; i<=degree; ++i) {
      if (std::abs(terms[i]) < std::numeric_limits<double>::epsilon()) continue;
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
    if (std::abs(terms[degree] - unity) == std::numeric_limits<double>::epsilon()) monic = true;
    if (std::abs(terms[0]) < std::numeric_limits<double>::epsilon()) homogeneous = true;
  }

  template<> 
  bool Polynomial<Rational>::is_null() const
  {
    if (terms.empty()) return true;
    unsigned int i;
    for(i=0; i<=degree; ++i) {
      if (!terms[i].is_null()) return false;
    }
    return true;
  }

  template<class kind> 
  bool Polynomial<kind>::is_null() const
  {
    if (terms.empty()) return true;
    unsigned int i;
    for(i=0; i<=degree; ++i) {
      if (std::abs(terms[i]) > std::numeric_limits<double>::epsilon()) return false;
    }
    return true;
  }

  template<class kind>
  std::ostream& operator <<(std::ostream& s,const Polynomial<kind>& source)
  {
    unsigned int i;

    if (source.terms[source.degree] > 0) {
      if (source.monic) {
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
    if (!source.homogeneous) {
      if (source.terms[0] > 0) {
        s << "+ " << source.terms[0];
      }
      else {
        s << "- " << -source.terms[0]; 		
      }
    } 	
    return s;
  }

  template<>
  bool operator ==(const Polynomial<Rational>& p1,const Polynomial<Rational>& p2)
  {
    if (p1.degree != p2.degree) return false;
    unsigned int i;
    Rational q;
    for(i=0; i<p1.degree; ++i) {
      q = p1.terms[i] - p2.terms[i];
      if (!q.is_null()) return false;
    }
    return true;
  }

  template<class kind>
  bool operator ==(const Polynomial<kind>& p1,const Polynomial<kind>& p2)
  {
    if (p1.degree != p2.degree) return false;
    unsigned int i;
    for(i=0; i<p1.degree; ++i) {
      if (std::abs(p1.terms[i] - p2.terms[i]) > std::numeric_limits<double>::epsilon()) return false;
    }
    return true;
  }

  template<class kind>
  Polynomial<kind> operator -(const Polynomial<kind>& p1,const Polynomial<kind>& p2)
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
    Polynomial<kind> output(new_terms);
    output.simplify();
    return output;
  }

  template<class kind>
  Polynomial<kind> operator +(const Polynomial<kind>& p1,const Polynomial<kind>& p2)
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
    Polynomial<kind> output(new_terms);
    output.simplify();
    return output;
  }

  template<class kind>
  Polynomial<kind> operator *(kind alpha,const Polynomial<kind>& p)
  {
    unsigned int i;
    std::vector<kind> new_terms;

    for(i=0; i<=p.degree; ++i) {
      new_terms.push_back(alpha*p.terms[i]);
    }
    Polynomial<kind> output(new_terms);
    output.simplify();
    return output;
  }

  template<class kind>
  Polynomial<kind> operator *(const Polynomial<kind>& p1,const Polynomial<kind>& p2)  
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
    Polynomial<kind> output(new_terms);
    output.simplify();
    return output;
  } 
} 
#endif





