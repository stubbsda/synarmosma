OBJECTS = global.o random.o cell.o nexus.o schema.o group.o graph.o geometry.o word.o rational.o\
variety_wrapper.o polynomial_wrapper.o proposition.o propositional_system.o logic_graph.o edge.o\
vertex.o multitime.o binary_matrix.o matrix_wrapper.o directed_graph.o functional_equation_wrapper.o\
homology.o homotopy.o poset.o lattice.o complex_graph.o

CXX_FLAGS += $(OPENMP) -fPIC #-DDISCRETE

#CXX_FLAGS += $(CXX_OPT)
CXX_FLAGS += $(DEBUG) 

#LD_FLAGS += $(CXX_OPT)
LD_FLAGS += $(DEBUG)

LD_FLAGS += $(OPENMP) -shared

UNAME = $(shell uname)

ifeq ($(UNAME),Darwin)
  LD_FLAGS += -install_name $(BIBLIOTHEK)/lib/libsynarmosma.so
endif 

LIBS = $(LAPACK) $(BOOST) -lntl -lm 

install: synarmosma
	mkdir -p $(BIBLIOTHEK)/lib
	install -p libsynarmosma.so $(BIBLIOTHEK)/lib/
	mkdir -p $(BIBLIOTHEK)/include/synarmosma
	install -p -m 444 *.h $(BIBLIOTHEK)/include/synarmosma/

synarmosma: $(OBJECTS) 
	$(CXX) $(LD_FLAGS) -o libsynarmosma.so $(OBJECTS) $(LIBS)   

schema.o: schema.cpp schema.h global.h
	$(CXX) $(CXX_FLAGS) -c schema.cpp

lattice.o: lattice.cpp lattice.h poset.h global.h
	$(CXX) $(CXX_FLAGS) -c lattice.cpp

poset.o: poset.cpp poset.h global.h
	$(CXX) $(CXX_FLAGS) -c poset.cpp

directed_graph.o: directed_graph.cpp directed_graph.h schema.h edge.h global.h
	$(CXX) $(CXX_FLAGS) -c directed_graph.cpp

rational.o: global.h rational.h rational.cpp
	$(CXX) $(CXX_FLAGS) -c rational.cpp

matrix_wrapper.o: matrix_wrapper.cpp matrix.cpp matrix.h global.h
	$(CXX) $(CXX_FLAGS) -c matrix_wrapper.cpp

proposition.o: proposition.cpp proposition.h global.h
	$(CXX) $(CXX_FLAGS) -c proposition.cpp

word.o: word.cpp word.h global.h
	$(CXX) $(CXX_FLAGS) -c word.cpp

homology.o: homology.h homology.cpp nexus.h matrix.h group.h global.h
	$(CXX) $(CXX_FLAGS) -c homology.cpp

homotopy.o: homotopy.h homotopy.cpp nexus.h matrix.h group.h global.h
	$(CXX) $(CXX_FLAGS) -c homotopy.cpp

binary_matrix.o: binary_matrix.cpp binary_matrix.h global.h
	$(CXX) $(CXX_FLAGS) -c binary_matrix.cpp

edge.o: global.h edge.h edge.cpp
	$(CXX) $(CXX_FLAGS) -c edge.cpp

complex_graph.o: global.h complex_graph.h complex_graph.cpp
	$(CXX) $(CXX_FLAGS) -c complex_graph.cpp

logic_graph.o: global.h graph.h logic_graph.h proposition.h propositional_system.h logic_graph.cpp
	$(CXX) $(CXX_FLAGS) -c logic_graph.cpp

polynomial_wrapper.o: global.h rational.h polynomial.h polynomial.cpp polynomial_wrapper.cpp
	$(CXX) $(CXX_FLAGS) -c polynomial_wrapper.cpp

variety_wrapper.o: global.h rational.h variety.h variety.cpp variety_wrapper.cpp
	$(CXX) $(CXX_FLAGS) -c variety_wrapper.cpp

functional_equation_wrapper.o: global.h rational.h variety.h polynomial.h functional_equation.h functional_equation.cpp functional_equation_wrapper.cpp
	$(CXX) $(CXX_FLAGS) -c functional_equation_wrapper.cpp

propositional_system.o: global.h graph.h proposition.h propositional_system.h propositional_system.cpp
	$(CXX) $(CXX_FLAGS) -c propositional_system.cpp

multitime.o: global.h multitime.h multitime.cpp
	$(CXX) $(CXX_FLAGS) -c multitime.cpp

vertex.o: global.h proposition.h vertex.h vertex.cpp
	$(CXX) $(CXX_FLAGS) -c vertex.cpp

graph.o: graph.cpp graph.h edge.h schema.h complex_graph.h matrix.h binary_matrix.h global.h
	$(CXX) $(CXX_FLAGS) -c graph.cpp

geometry.o: geometry.cpp geometry.h global.h 
	$(CXX) $(CXX_FLAGS) -c geometry.cpp

cell.o: cell.cpp cell.h global.h
	$(CXX) $(CXX_FLAGS) -c cell.cpp

nexus.o: nexus.cpp nexus.h cell.h schema.h group.h matrix.h word.h global.h
	$(CXX) $(CXX_FLAGS) -c nexus.cpp

group.o: group.cpp group.h word.h global.h
	$(CXX) $(CXX_FLAGS) -c group.cpp

global.o: global.cpp global.h
	$(CXX) $(CXX_FLAGS) -c global.cpp

random.o: random.cpp global.h
	$(CXX) $(CXX_FLAGS) -c random.cpp

clean:
	rm -f $(OBJECTS)
	rm -f *~
	rm -f libsynarmosma.so
	rm -f $(BIBLIOTHEK)/lib/libsynarmosma.so
	rm -rf $(BIBLIOTHEK)/include/synarmosma











