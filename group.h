#include "word.h"

#ifndef _grouph
#define _grouph

namespace SYNARMOSMA {
  // A class for a combinatorial group presentation
  class Group {
   private:
    unsigned int ngenerator;
    std::vector<Word> relations;
    bool abelian;
    bool finite;
    bool solvable;
    bool free;
    unsigned int cardinality;
    unsigned int rank;
    std::vector<unsigned int> torsion;

    void allocate(unsigned int);
    void compute_rank();
    Group abelianize() const;
   public:
    Group();
    Group(int);
    Group(int,const std::vector<Word>&);
    Group(int,int);
    Group(unsigned int,const std::vector<unsigned int>&);
    Group(const Group&);
    Group(const std::string&,int);
    ~Group();
    Group& operator =(const Group&);
    inline int get_rank() const {return rank;};
    unsigned int implied_generators() const;
    void initialize(unsigned int,const std::vector<Word>&);
    void initialize(unsigned int,const std::vector<unsigned int>&);
    void reduce();
    void clear();
    void create_random();
    bool consistent() const;
    void serialize(std::ofstream&) const;
    void deserialize(std::ifstream&);
    std::string compact_form() const;
    friend std::ostream& operator <<(std::ostream&,const Group&);
    friend class Homotopy;
  };
}
#endif
