# Which C++ compiler to use, default is Gnu g++
CXX         = g++
# General flags for C++ compilation, assumes compiler with OpenMP support
CXX_FLAGS   = -std=c++14 -Wall -march=native -fPIC -fopenmp
# For debug mode
DEBUG       = -g -DDEBUG
# Optimization values for the g++ on an x86-64 platform
OPT         = -O3 -fstrict-aliasing -ftree-vectorize -funroll-loops
# Choose either the debug or optimized version of the library
#CXX_FLAGS  += $(OPT)
CXX_FLAGS  += $(DEBUG)
# Use the 64 bit integer versions of the vertex energy and spatial geometry?
#CXX_FLAGS += -DDISCRETE
# Standard linking flags, shouldn't need to be changed
LD_FLAGS    = -Wall -shared -fopenmp
# Debug or optimized version for linking...
#LD_FLAGS   += $(OPT)
LD_FLAGS   += $(DEBUG)
# Name of the LAPACK and BLAS libraries on this platform
LAPACK      = -llapack -lblas
# The Boost library on this platform
BOOST       = -lboost_system
# All the necessary libraries: LAPACK, Boost, NTL and GMP 
LIBS        = $(LAPACK) $(BOOST) -lntl -lgmp -lm
# Where to install the Synarmosma library, requires root privileges
INSTALL_DIR = $(HOME)/fabrica/local
# End of user-modifiable parameters!

OBJECTS = global.o random.o cell.o nexus.o schema.o group.o graph.o geometry.o word.o rational.o\
variety_wrapper.o polynomial_wrapper.o integer_polynomial_wrapper.o proposition.o propositional_system.o\
vertex.o multitime.o binary_matrix.o integer_matrix_wrapper.o matrix_wrapper.o logic_graph.o edge.o\
functional_equation_wrapper.o homology.o homotopy.o poset.o lattice.o pseudograph.o solver_wrapper.o\
directed_graph.o

UNAME = $(shell uname)

ifeq ($(UNAME),Darwin)
  RPATH += -install_name $(INSTALL_DIR)/lib/libsynarmosma.so
endif 

ifeq ($(UNAME),Linux)
  RPATH = -Wl,-rpath $(INSTALL_DIR)/lib
endif 

install: synarmosma
	mkdir -p $(INSTALL_DIR)/lib
	install -p libsynarmosma.so $(INSTALL_DIR)/lib/
	mkdir -p $(INSTALL_DIR)/include/synarmosma
	install -p -m 444 *.h $(INSTALL_DIR)/include/synarmosma/

test: install
	$(CXX) -std=c++14 -Wall -march=native $(RPATH) -I$(INSTALL_DIR)/include -L$(INSTALL_DIR)/lib -o unit_test unit_testing.cpp -lsynarmosma $(LIBS)
	./unit_test

synarmosma: $(OBJECTS) 
	$(CXX) $(RPATH) $(LD_FLAGS) -o libsynarmosma.so $(OBJECTS) $(LIBS)   

schema.o: schema.cpp schema.h global.h
	$(CXX) $(CXX_FLAGS) -c schema.cpp

lattice.o: lattice.cpp lattice.h poset.h random.h global.h
	$(CXX) $(CXX_FLAGS) -c lattice.cpp

poset.o: poset.cpp poset.h random.h global.h
	$(CXX) $(CXX_FLAGS) -c poset.cpp

directed_graph.o: directed_graph.cpp directed_graph.h schema.h pseudograph.h matrix.h binary_matrix.h edge.h random.h global.h
	$(CXX) $(CXX_FLAGS) -c directed_graph.cpp

rational.o: global.h rational.h rational.cpp
	$(CXX) $(CXX_FLAGS) -c rational.cpp

solver_wrapper.o: solver_wrapper.cpp solver.cpp solver.h matrix.h random.h global.h
	$(CXX) $(CXX_FLAGS) -c solver_wrapper.cpp

integer_matrix_wrapper.o: integer_matrix_wrapper.cpp integer_matrix.cpp integer_matrix.h random.h global.h
	$(CXX) $(CXX_FLAGS) -c integer_matrix_wrapper.cpp

matrix_wrapper.o: matrix_wrapper.cpp matrix.cpp matrix.h random.h global.h
	$(CXX) $(CXX_FLAGS) -c matrix_wrapper.cpp

proposition.o: proposition.cpp proposition.h random.h global.h
	$(CXX) $(CXX_FLAGS) -c proposition.cpp

word.o: word.cpp word.h random.h global.h
	$(CXX) $(CXX_FLAGS) -c word.cpp

homology.o: homology.h homology.cpp nexus.h integer_matrix.h binary_matrix.h graph.h cell.h schema.h random.h global.h
	$(CXX) $(CXX_FLAGS) -c homology.cpp

homotopy.o: homotopy.h homotopy.cpp nexus.h cell.h schema.h graph.h edge.h group.h word.h random.h global.h
	$(CXX) $(CXX_FLAGS) -c homotopy.cpp

binary_matrix.o: binary_matrix.cpp binary_matrix.h random.h global.h
	$(CXX) $(CXX_FLAGS) -c binary_matrix.cpp

edge.o: global.h edge.h edge.cpp
	$(CXX) $(CXX_FLAGS) -c edge.cpp

pseudograph.o: global.h pseudograph.h pseudograph.cpp
	$(CXX) $(CXX_FLAGS) -c pseudograph.cpp

logic_graph.o: global.h random.h edge.h schema.h pseudograph.h matrix.h binary_matrix.h graph.h logic_graph.h proposition.h propositional_system.h logic_graph.cpp
	$(CXX) $(CXX_FLAGS) -c logic_graph.cpp

polynomial_wrapper.o: global.h rational.h polynomial.h random.h polynomial.cpp polynomial_wrapper.cpp
	$(CXX) $(CXX_FLAGS) -c polynomial_wrapper.cpp

integer_polynomial_wrapper.o: global.h integer_polynomial.h random.h integer_polynomial.cpp integer_polynomial_wrapper.cpp
	$(CXX) $(CXX_FLAGS) -c integer_polynomial_wrapper.cpp

variety_wrapper.o: global.h rational.h variety.h random.h variety.cpp variety_wrapper.cpp
	$(CXX) $(CXX_FLAGS) -c variety_wrapper.cpp

functional_equation_wrapper.o: global.h rational.h variety.h polynomial.h functional_equation.h functional_equation.cpp functional_equation_wrapper.cpp
	$(CXX) $(CXX_FLAGS) -c functional_equation_wrapper.cpp

propositional_system.o: global.h graph.h proposition.h propositional_system.h propositional_system.cpp binary_matrix.h directed_graph.h
	$(CXX) $(CXX_FLAGS) -c propositional_system.cpp

multitime.o: global.h multitime.h multitime.cpp
	$(CXX) $(CXX_FLAGS) -c multitime.cpp

vertex.o: global.h vertex.h random.h vertex.cpp
	$(CXX) $(CXX_FLAGS) -c vertex.cpp

graph.o: graph.cpp graph.h edge.h schema.h pseudograph.h matrix.h binary_matrix.h random.h global.h
	$(CXX) $(CXX_FLAGS) -c graph.cpp

geometry.o: geometry.cpp geometry.h matrix.h random.h global.h 
	$(CXX) $(CXX_FLAGS) -c geometry.cpp

cell.o: cell.cpp cell.h global.h
	$(CXX) $(CXX_FLAGS) -c cell.cpp

nexus.o: nexus.cpp nexus.h cell.h schema.h graph.h edge.h random.h matrix.h binary_matrix.h pseudograph.h global.h
	$(CXX) $(CXX_FLAGS) -c nexus.cpp

group.o: group.cpp group.h word.h random.h global.h
	$(CXX) $(CXX_FLAGS) -c group.cpp

global.o: global.cpp random.h global.h
	$(CXX) $(CXX_FLAGS) -c global.cpp

random.o: random.cpp random.h global.h
	$(CXX) $(CXX_FLAGS) -c random.cpp

clean:
	rm -f $(OBJECTS)
	rm -f *~
	rm -f unit_test
	rm -f libsynarmosma.so
	rm -f $(BIBLIOTHEK)/lib/libsynarmosma.so
	rm -rf $(BIBLIOTHEK)/include/synarmosma











