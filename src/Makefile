include ../Makefile.config

OBJECTS = global.o random.o cell.o nexus.o schema.o group.o graph.o geometry_wrapper.o word.o rational.o\
variety_wrapper.o polynomial_wrapper.o integer_polynomial_wrapper.o proposition.o propositional_system.o\
vertex_wrapper.o multitime_wrapper.o binary_matrix.o integer_matrix_wrapper.o matrix_wrapper.o logic_graph.o\
functional_equation.o homology.o homotopy.o markov_chain.o lattice.o pseudograph.o solver_wrapper.o edge.o\
directed_graph.o poset.o

LD_FLAGS += -Wl,-rpath $(INSTALL_DIR)/lib

synarmosma: $(OBJECTS)
	$(CXX) $(LD_FLAGS) -o libsynarmosma.so $(OBJECTS) $(LIBS)   

schema.o: schema.cxx schema.h global.h
	$(CXX) $(CXX_FLAGS) -c schema.cxx

lattice.o: lattice.cxx lattice.h poset.h random.h global.h
	$(CXX) $(CXX_FLAGS) -c lattice.cxx

poset.o: poset.cxx poset.h graph.h random.h global.h
	$(CXX) $(CXX_FLAGS) -c poset.cxx

markov_chain.o: markov_chain.cxx markov_chain.h matrix.h random.h global.h
	$(CXX) $(CXX_FLAGS) -c markov_chain.cxx

directed_graph.o: directed_graph.cxx directed_graph.h graph.h random.h global.h
	$(CXX) $(CXX_FLAGS) -c directed_graph.cxx

rational.o: rational.cxx rational.h global.h
	$(CXX) $(CXX_FLAGS) -c rational.cxx

solver_wrapper.o: solver_wrapper.cxx solver.cxx solver.h graph.h matrix.h random.h global.h
	$(CXX) $(CXX_FLAGS) -c solver_wrapper.cxx

integer_matrix_wrapper.o: integer_matrix_wrapper.cxx integer_matrix.cxx integer_matrix.h global.h
	$(CXX) $(CXX_FLAGS) -c integer_matrix_wrapper.cxx

matrix_wrapper.o: matrix_wrapper.cxx matrix.cxx matrix.h global.h
	$(CXX) $(CXX_FLAGS) -c matrix_wrapper.cxx

proposition.o: proposition.cxx proposition.h random.h global.h
	$(CXX) $(CXX_FLAGS) -c proposition.cxx

word.o: word.cxx word.h random.h global.h
	$(CXX) $(CXX_FLAGS) -c word.cxx

homology.o: homology.cxx homology.h integer_matrix.h binary_matrix.h nexus.h global.h
	$(CXX) $(CXX_FLAGS) -c homology.cxx

homotopy.o: homotopy.cxx homotopy.h group.h word.h nexus.h random.h global.h
	$(CXX) $(CXX_FLAGS) -c homotopy.cxx

binary_matrix.o: binary_matrix.cxx binary_matrix.h random.h global.h
	$(CXX) $(CXX_FLAGS) -c binary_matrix.cxx

edge.o: edge.cxx edge.h global.h  
	$(CXX) $(CXX_FLAGS) -c edge.cxx

pseudograph.o: pseudograph.cxx pseudograph.h edge.h global.h  
	$(CXX) $(CXX_FLAGS) -c pseudograph.cxx

logic_graph.o: logic_graph.cxx logic_graph.h propositional_system.h graph.h edge.h random.h global.h 
	$(CXX) $(CXX_FLAGS) -c logic_graph.cxx

polynomial_wrapper.o: polynomial_wrapper.cxx polynomial.cxx polynomial.h rational.h random.h global.h
	$(CXX) $(CXX_FLAGS) -c polynomial_wrapper.cxx

integer_polynomial_wrapper.o: integer_polynomial_wrapper.cxx integer_polynomial.cxx integer_polynomial.h group.h rational.h global.h    
	$(CXX) $(CXX_FLAGS) -c integer_polynomial_wrapper.cxx

variety_wrapper.o: variety_wrapper.cxx variety.cxx variety.h graph.h rational.h random.h global.h
	$(CXX) $(CXX_FLAGS) -c variety_wrapper.cxx

functional_equation.o: functional_equation.cxx functional_equation.h variety.h polynomial.h rational.h global.h
	$(CXX) $(CXX_FLAGS) -c functional_equation.cxx

propositional_system.o: propositional_system.cxx propositional_system.h proposition.h directed_graph.h binary_matrix.h global.h
	$(CXX) $(CXX_FLAGS) -c propositional_system.cxx

multitime_wrapper.o: multitime_wrapper.cxx multitime.cxx multitime.h global.h
	$(CXX) $(CXX_FLAGS) -c multitime_wrapper.cxx

vertex_wrapper.o: vertex_wrapper.cxx vertex.cxx vertex.h global.h   
	$(CXX) $(CXX_FLAGS) -c vertex_wrapper.cxx

graph.o: graph.cxx graph.h edge.h schema.h pseudograph.h matrix.h binary_matrix.h integer_polynomial.h random.h global.h
	$(CXX) $(CXX_FLAGS) -c graph.cxx

geometry_wrapper.o: geometry_wrapper.cxx geometry.cxx geometry.h matrix.h random.h global.h 
	$(CXX) $(CXX_FLAGS) -c geometry_wrapper.cxx

cell.o: cell.cxx cell.h global.h
	$(CXX) $(CXX_FLAGS) -c cell.cxx

nexus.o: nexus.cxx nexus.h graph.h cell.h schema.h global.h
	$(CXX) $(CXX_FLAGS) -c nexus.cxx

group.o: group.cxx group.h word.h random.h global.h
	$(CXX) $(CXX_FLAGS) -c group.cxx

global.o: global.cxx global.h
	$(CXX) $(CXX_FLAGS) -c global.cxx

random.o: random.cxx random.h global.h
	$(CXX) $(CXX_FLAGS) -c random.cxx

clean:
	rm -f $(OBJECTS)
	rm -f *~
	rm -f libsynarmosma.so

