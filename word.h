#include "global.h"

#ifndef _wordh
#define _wordh

class Word {
 private:
  std::vector<std::pair<unsigned int,int> > content;
  // The number of letters in the alphabet for words in a group presentation
  unsigned int NL;

  void initialize(unsigned int);
  void initialize(const std::vector<unsigned int>&,const std::vector<int>&);
  void clear();
 public:
  Word(unsigned int);
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
  Word reduce(int,const std::set<unsigned int>&,const unsigned int*) const;
  unsigned int size() const;
  void permute(unsigned int,Word&) const;
  bool trivial() const;
  bool alias() const;
  bool legal() const;
  bool empty() const;
  void serialize(std::ofstream&) const;
  void deserialize(std::ifstream&);
  friend bool operator ==(const Word&,const Word&);
  friend bool operator !=(const Word&,const Word&);
  friend Word operator *(const Word&,const Word&);
  friend std::ostream& operator <<(std::ostream&,const Word&);
  friend class Group;
  friend class Nexus;
  friend class Spacetime;
};
#endif
