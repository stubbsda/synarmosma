#include "random.h"

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

    /// This method initializes the word to the form \f$w = x_0^{k_0}\cdots x_{n-1}^{k_{n-1}}\f$ where \f$n\f$ is the argument and the \f$k_i\f$ are randomly chosen between -10 and 10, excluding zero.
    void initialize(unsigned int);
    /// This method initializes the word to the form \f$w = x_n^m\f$ where \f$n\f$ is the first argument and \f$m\ne 0\f$ the second argument.
    void initialize(unsigned int,int);
    /// This method initializes the word based on the first argument being a vector of letter indices and the second argument the exponents for these same letters.
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
    void clear();
    /// This method appends a new "letter" to the end of the word.
    void append(const std::pair<unsigned int,int>&);
    /// This method computes the set of letters used in this word and returns the cardinality of this set.
    unsigned int get_alphabet(std::set<unsigned int>&) const;
    /// This operator has the same effect as calling the invert method on this instance of the Word class.
    Word operator !() const;
    /// This method returns the inverse of the word, i.e. a word with the order of letters reversed and the sign of all the exponents reversed.
    Word invert() const;
    /// This method mutates the word, i.e. a random letter is selected from the word and converted to another letter; if the word contains just one letter then the exponent is altered to a randomly chosen value between -10 and 10, excluding zero.
    Word mutate() const;
    /// This method eliminates redundant elements of the word, fusing together adjacent letters which share the same index and dropping letters whose exponent is zero.
    Word normalize() const;
    /// This method swaps every occurrence in the word of the second index by the first index, reversing the exponent's sign when the final argument is true, and reducing the index values by one when they are greater than the second argument. 
    Word swap(unsigned int,unsigned int,bool = false) const;
    /// This method creates a new word from an existing one in a similar manner to the Word::normalize method, identifying pairs of adjacent letters with the same index and opposite but equal exponents, so that the pair can be eliminated.
    Word reduce() const;
    /// This method creates a new word from the existing one by eliminating all of the letters whose index lies in the first argument, using the array of offset indices in the second argument. 
    Word reduce(const std::set<unsigned int>&,const unsigned int*) const;
    /// This method returns the length of the vector "content".
    unsigned int length() const;
    /// This method returns true when the "content" vector is empty.
    bool empty() const;
    /// This method carries out a cyclic permutation based on its unique argument \f$k\f$, so that (ignoring exponents) if \f$w = x_0\cdots x_{n-1}\f$ then \f$w' = x_k x_{k+1}\cdots x_{n-1} x_0 x_1 \cdots x_{k-1}\f$. 
    Word permute(unsigned int) const;
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
  };

  inline void Word::clear()
  {
    content.clear();
  }

  inline void Word::append(const std::pair<unsigned int,int>& term)
  {
    content.push_back(term);
  }

  inline unsigned int Word::length() const
  {
    return content.size();
  }

  inline bool Word::empty() const
  {
    return content.empty();
  }
}
#endif
