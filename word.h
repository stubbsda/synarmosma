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
   public:
    /// The default constructor that leaves the instance empty.
    Word();
    /// This constructor creates a word of the form \f$w = x_0^{k_0}\cdots x_{n-1}^{k_{n-1}}\f$ where \f$n\f$ is the argument and the \f$k_i\f$ are randomly chosen between -10 and 10, excluding zero.
    Word(unsigned int);
    /// This constructor builds a word with a single element, \f$w = x_n^m\f$ where \f$n\f$ is the first argument and \f$m\ne 0\f$ the second argument.
    Word(unsigned int,int);
    /// This constructor accepts a string of lower-case and upper-case letters (representing negative exponents) to build a word; thus an expression like "aaBcAf" would become \f$x_0^2 x_1^{-1}x_2 x_0^{-1}x_5\f$.
    Word(const std::string&);
    /// The standard copy constructor that copies over the value of the "content" property.
    Word(const Word&);
    /// The overloaded assignment operator that copies over the value of the "content" property.
    Word& operator =(const Word&);
    /// The destructor for this class, which does nothing.
    ~Word();
    /// This method calls clear on the vector "content" and so empties the word of its contents.
    inline void clear() {content.clear();}; 
    /// This method computes the set of letters used in this word and returns the cardinality of this set.
    unsigned int get_alphabet(std::set<unsigned int>&) const;
    /// This operator has the same effect as calling the invert method on this instance of the Word class.
    Word operator !() const;
    /// This method returns the inverse of the word, i.e. a word with the order of letters reversed and the sign of all the exponents reversed.
    Word invert() const;
    Word mutate() const;
    Word normalize() const;
    Word swap(unsigned int,unsigned int,bool = false) const;
    void free_reduce();
    Word reduce(const std::set<unsigned int>&,const unsigned int*) const;
    /// This method returns the length of the vector "content".
    inline unsigned int length() const {return content.size();};
    /// This method returns true when the "content" vector is empty.
    inline bool empty() const {return content.empty();};
    void permute(unsigned int,Word&) const;
    /// This method determines if the word consists of a single letter whose exponent is \f$\pm 1\f$ and returns true if this is so.
    bool trivial() const;
    /// This method returns true when the word has the form \f$w = x_p^n x_q^m\f$ with \f$p \ne q\f$ and \f$|n| = |m| =1\f$.
    bool alias() const;
    /// This method verifies that none of the exponents in the "content" vector are equal to zero, returning true if this is so.
    bool legal() const;
    /// This method determines if the lowest index letter in the word occurs in at least two distinct locations with exponents whose sign differs, returning false in this case and true otherwise.
    bool homogeneous() const;
    /// This method writes the instance properties to a binary disk file and returns the number of bytes written to the file.
    int serialize(std::ofstream&) const;
    /// This method calls the clear() method on the instance and then reads the properties from a binary disk file and returns the number of bytes read.
    int deserialize(std::ifstream&);
    /// This overloaded operator tests if two words are equal, i.e. the each letter and its exponent is the same in the two words.
    friend bool operator ==(const Word&,const Word&);
    /// This overloaded operator simply uses the "==" operator to test the two words for equality and reverses the output.
    friend bool operator !=(const Word&,const Word&);
    /// This overloaded operator concatenates the two words, calls the Word::normalize method on this product word and returns the output.
    friend Word operator *(const Word&,const Word&);
    /// This overloading of the ostream operator will print the word in a nice human-readable format with the letters as "x" with a subscript.
    friend std::ostream& operator <<(std::ostream&,const Word&);
    /// This function accepts two words as argument and computes the intersection of their respective alphabet, returning the cardinality of this intersection.
    friend int affinity(const Word&,const Word&,std::set<unsigned int>&);
    friend class Group;
    friend class Homotopy;
  };
}
#endif
