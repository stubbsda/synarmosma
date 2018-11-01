#include "word.h"

#ifndef _grouph
#define _grouph

namespace SYNARMOSMA {
  // A class for a combinatorial group presentation
  class Group {
   protected:
    unsigned int ngenerator = 0;
    unsigned int cardinality = 0;
    unsigned int rank = 0;
    bool abelian = false;
    bool finite = false;
    bool solvable = false;
    bool free = false;
    bool braid = false;
    std::vector<unsigned int> torsion;
    std::vector<Word> relations;

    void allocate(unsigned int);
    void compute_rank();
    Group abelianize() const;
   public:
    Group();
    Group(unsigned int);
    Group(unsigned int,const std::vector<Word>&);
    Group(unsigned int,unsigned int);
    Group(unsigned int,const std::vector<unsigned int>&);
    Group(const Group&);
    Group(const std::string&,unsigned int);
    ~Group();
    Group& operator =(const Group&);
    inline unsigned int get_rank() const {return rank;};
    int implied_generators() const;
    void initialize(unsigned int,const std::vector<Word>&);
    void initialize(unsigned int,const std::vector<unsigned int>&);
    bool equivalent(const Word&,const Word&) const;
    void reduce();
    void clear();
    void create_random();
    bool consistent() const;
    int serialize(std::ofstream&) const;
    int deserialize(std::ifstream&);
    std::string compact_form() const;
    friend std::ostream& operator <<(std::ostream&,const Group&);
    friend class Homotopy;
  };
}
#endif
