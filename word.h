#include "global.h"

#ifndef _wordh
#define _wordh

namespace SYNARMOSMA {
  class Word {
   private:
    std::vector<std::pair<unsigned int,int> > content;
    // The number of letters in the alphabet for words in a group presentation
    unsigned int NL = 0;

    void initialize(unsigned int);
    void initialize(unsigned int,unsigned int,int);
    void initialize(const std::vector<unsigned int>&,const std::vector<int>&);
    inline void clear() {NL = 0; content.clear();};
   public:
    Word();
    Word(int);
    Word(int,const std::string&);
    Word(int,int);
    Word(int,int,int);
    Word(const Word&);
    Word& operator =(const Word&);
    ~Word();
    Word operator !() const;
    Word invert() const;
    Word mutate() const;
    Word normalize() const;
    Word swap(int,int,bool) const;
    void free_reduce();
    Word reduce(unsigned int,const std::set<unsigned int>&,const unsigned int*) const;
    inline unsigned int length() const {return content.size();};
    inline bool empty() const {return content.empty();};
    void permute(int,Word&) const;
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
    friend class Group;
    friend class Homotopy;
  };
}
#endif
