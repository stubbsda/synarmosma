/*
  Copyright 2014 Daniel Stubbs

  This file is part of Synarmosma.

  Synarmosma is free software: you can redistribute it and/or modify 
  it under the terms of the GNU General Public License as published by 
  the Free Software Foundation, either version 3 of the License, or 
  (at your option) any later version.

  Synarmosma is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with Synarmosma.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "variety.h"
#include "polynomial.h"

#ifndef _functionaleqnh
#define _functionaleqnh

namespace SYNARMOSMA {
  template<class kind>
  class Functional_Equation;

  template<class kind>
  std::ostream& operator <<(std::ostream& s,const Functional_Equation<kind>&);

  template<class kind>
  class Functional_Equation {
   private:
    // We need a triple from Boost!
    std::vector<boost::tuple<Polynomial<kind>,Polynomial<kind>,unsigned int> > terms;
    Polynomial<kind> remainder;
    bool linear;
    bool homogeneous;
   
    void initialize(unsigned int);
    void analyze_file(std::vector<std::string>&,std::vector<std::string>&,std::vector<std::string>&);
   public:
    Functional_Equation();
    Functional_Equation(unsigned int);
    Functional_Equation(const char*);
    Functional_Equation(const Functional_Equation&);
    Functional_Equation& operator =(const Functional_Equation&);
    ~Functional_Equation();
    void serialize(std::ofstream&) const;
    void deserialize(std::ifstream&);
    Variety<unsigned int> reduce(unsigned int);  
    friend std::ostream& operator << <>(std::ostream& s,const Functional_Equation<kind>&);
  };

  template<class kind>
  std::ostream& operator <<(std::ostream& s,const Functional_Equation<kind>& source)
  {
    unsigned int i;
    boost::tuple<Polynomial<kind>,Polynomial<kind>,unsigned int> trio;
    for(i=0; i<source.terms.size(); ++i) {
      trio = source.terms[i];
      s << "(" << boost::get<0>(trio) << ")*F(" << boost::get<1>(trio) << ")^" << boost::get<2>(trio) << " +" << std::endl; 
    }
    s << source.remainder << " = 0";
    return s;
  }
}
#endif

