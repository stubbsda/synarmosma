#include "matrix.h"
#include "random.h"

namespace SYNARMOSMA {
  /// A class representing a general discrete-time Markov chain for a system with a finite number of distinct states.
  class Markov_Chain {
   private:
     /// This property stores the dimension of the Markov chain, i.e. the number of distinct states in the system, 
     /// which must be at least two. Once it is set in the constructor it cannot be changed.
     unsigned int N = 0;
     /// This property stores the Markov chain's transition matrix, a square matrix of dimension two or more, all of 
     /// whose entries are non-negative and such that the sum of each row is one. 
     Matrix<double>* transition_matrix;

     /// This method multiplies the state vector of probabilities which is its first argument by the transition matrix and assigns the second argument to this product.  
     void multiply(const std::vector<double>&,std::vector<double>&) const;
   public:
     /// The constructor for this class, whose two arguments are the values of Markov_Chain::N and the extent of the initial mutation for the transition matrix. After the Markov_Chain::transition_matrix property is allocated and initialized as the identity matrix, the constructor calls Markov_Chain::mutate with the value of the second argument.
     Markov_Chain(unsigned int,double = 0.0);
     /// The destructor for this class which restores the memory for the Markov_Chain::transition_matrix property, assuming it has been allocated.
     ~Markov_Chain();
     /// This method changes the value of the transition matrix while still respecting its constraints; the severity of the mutation is determined by the method's positive argument (if it is zero, the method exits immediately) - the greater the argument the more pronounced the mutation. After the changes are made the Markov_Chain::consistent method is called to verify coherence of the transition matrix. 
     void mutate(double);
     /// This method sets a given row of the transition matrix, specified by the second argument, to the vector of probabilities that is the method's first argument. After the changes are made the Markov_Chain::consistent method is called to verify coherence of the transition matrix. 
     void set_transition(const std::vector<double>&,unsigned int);
     /// This method ensures that the properties of this instance of the class are coherent: the dimension is greater than one, the transition matrix has the correct number of rows and columns and all of its elements satisfy the constraints for a transition matrix (i.e. all elements are non-negative and each row's elements sum to unity). If this is so it returns true, false otherwise. 
     bool consistent() const;
     /// This method computes the next or later state, depending on the third argument, as a vector of probabilities assigned to the second argument, given that it is in the state specified by the method's first argument; it calls the Markov_Chain::multiply method to accomplish this.
     void get_state(const std::vector<double>&,std::vector<double>&,int = 1) const;
     /// This method returns the next or later state, depending on the second argument, for the system by a random number, given that it is in the state specified by the method's first argument; it calls the Markov_Chain::multiply method to accomplish this.
     int get_state(const std::vector<double>&,int = 1) const;
     /// This method returns the next state for the system by a random number, given that it is in the state specified by the method's argument.
     int get_state(unsigned int) const;
  };
}