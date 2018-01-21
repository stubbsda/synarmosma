#include "global.h"

#ifndef _wordh
#define _wordh

namespace SYNARMOSMA {
  class Word {
   private:
    std::vector<std::pair<unsigned int,int> > content;
    // The number of letters in the alphabet for words in a group presentation
    unsigned int NL;

    void initialize(unsigned int);
    void initialize(const std::vector<unsigned int>&,const std::vector<int>&);
    inline void clear() {content.clear();};
   public:
    Word();
    Word(unsigned int);
    Word(unsigned int,const std::string&);
    Word(unsigned int,unsigned int);
    Word(unsigned int,unsigned int,int);
    Word(const Word&);
    Word& operator =(const Word&);
    ~Word();
    Word operator !() const;
    Word invert() const;
    Word mutate() const;
    Word normalize() const;
    Word swap(unsigned int,unsigned int,bool) const;
    void free_reduce();
    Word reduce(int,const std::set<unsigned int>&,const unsigned int*) const;
    inline unsigned int length() const {return content.size();};
    inline bool empty() const {return content.empty();};
    void permute(unsigned int,Word&) const;
    bool trivial() const;
    bool alias() const;
    bool legal() const;
    bool homogeneous() const;
    int serialize(std::ofstream&) const;
    int deserialize(std::ifstream&);
    void initialize(unsigned int,unsigned int,int);
    void write2screen() const;
    friend bool operator ==(const Word&,const Word&);
    friend bool operator !=(const Word&,const Word&);
    friend Word operator *(const Word&,const Word&);
    friend std::ostream& operator <<(std::ostream&,const Word&);
    friend class Group;
    friend class Homotopy;
  };
}
#endif
