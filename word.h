#include "global.h"

#ifndef _wordh
#define _wordh

namespace SYNARMOSMA {
  /// A class representing a word in an alphabet, that is an expression involving the finite product of atomic letters, raised to integer powers. 
  class Word {
   private:
    /// The word is stored as a vector of pairs - the first element of the pair 
    /// is the letter, an unsigned integer, while the second element of the pair 
    /// is the letter's exponent, a non-zero integer. 
    std::vector<std::pair<unsigned int,int> > content;

    void initialize(unsigned int);
    void initialize(unsigned int,int);
    void initialize(const std::vector<unsigned int>&,const std::vector<int>&);
    inline void clear() {content.clear();}; 
   public:
    Word();
    Word(unsigned int);
    Word(unsigned int,int);
    Word(const std::string&);
    Word(const Word&);
    Word& operator =(const Word&);
    ~Word();
    unsigned int get_alphabet(std::set<unsigned int>&) const;
    Word operator !() const;
    Word invert() const;
    Word mutate() const;
    Word normalize() const;
    Word swap(unsigned int,unsigned int,bool = false) const;
    void free_reduce();
    Word reduce(unsigned int,const std::set<unsigned int>&,const unsigned int*) const;
    inline unsigned int length() const {return content.size();};
    inline bool empty() const {return content.empty();};
    void permute(unsigned int,Word&) const;
    bool trivial() const;
    bool alias() const;
    bool legal() const;
    bool homogeneous() const;
    int serialize(std::ofstream&) const;
    int deserialize(std::ifstream&);
    friend bool operator ==(const Word&,const Word&);
    friend bool operator !=(const Word&,const Word&);
    friend Word operator *(const Word&,const Word&);
    friend std::ostream& operator <<(std::ostream&,const Word&);
    friend int affinity(const Word&,const Word&,std::set<unsigned int>&);
    friend class Group;
    friend class Homotopy;
  };
}
#endif
