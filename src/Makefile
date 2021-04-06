include ../Makefile.config

OBJECTS = global.o random.o cell.o nexus.o schema.o group.o graph.o geometry_wrapper.o word.o rational.o\
variety_wrapper.o polynomial_wrapper.o integer_polynomial_wrapper.o proposition.o propositional_system.o\
vertex_wrapper.o multitime.o binary_matrix.o integer_matrix_wrapper.o matrix_wrapper.o logic_graph.o edge.o\
functional_equation_wrapper.o homology.o homotopy.o markov_chain.o lattice.o pseudograph.o solver_wrapper.o\
directed_graph.o poset.o

LD_FLAGS += -Wl,-rpath $(INSTALL_DIR)/lib

synarmosma: $(OBJECTS)
	$(CXX) $(LD_FLAGS) -o libsynarmosma.so $(OBJECTS) $(LIBS)   

schema.o: schema.cxx schema.h global.h
	$(CXX) $(CXX_FLAGS) -c schema.cxx

lattice.o: lattice.cxx lattice.h poset.h random.h global.h
	$(CXX) $(CXX_FLAGS) -c lattice.cxx

poset.o: poset.cxx poset.h random.h global.h
	$(CXX) $(CXX_FLAGS) -c poset.cxx

markov_chain.o: markov_chain.cxx markov_chain.h matrix.h random.h global.h
	$(CXX) $(CXX_FLAGS) -c markov_chain.cxx

directed_graph.o: directed_graph.cxx directed_graph.h schema.h pseudograph.h matrix.h binary_matrix.h edge.h random.h global.h
	$(CXX) $(CXX_FLAGS) -c directed_graph.cxx

rational.o: global.h rational.h rational.cxx
	$(CXX) $(CXX_FLAGS) -c rational.cxx

solver_wrapper.o: solver_wrapper.cxx solver.cxx solver.h matrix.h random.h global.h
	$(CXX) $(CXX_FLAGS) -c solver_wrapper.cxx

integer_matrix_wrapper.o: integer_matrix_wrapper.cxx integer_matrix.cxx integer_matrix.h random.h global.h
	$(CXX) $(CXX_FLAGS) -c integer_matrix_wrapper.cxx

matrix_wrapper.o: matrix_wrapper.cxx matrix.cxx matrix.h random.h global.h
	$(CXX) $(CXX_FLAGS) -c matrix_wrapper.cxx

proposition.o: proposition.cxx proposition.h random.h global.h
	$(CXX) $(CXX_FLAGS) -c proposition.cxx

word.o: word.cxx word.h random.h global.h
	$(CXX) $(CXX_FLAGS) -c word.cxx

homology.o: homology.h homology.cxx nexus.h integer_matrix.h binary_matrix.h graph.h cell.h schema.h random.h global.h
	$(CXX) $(CXX_FLAGS) -c homology.cxx

homotopy.o: homotopy.h homotopy.cxx nexus.h cell.h schema.h graph.h edge.h group.h word.h random.h global.h
	$(CXX) $(CXX_FLAGS) -c homotopy.cxx

binary_matrix.o: binary_matrix.cxx binary_matrix.h random.h global.h
	$(CXX) $(CXX_FLAGS) -c binary_matrix.cxx

edge.o: global.h edge.h edge.cxx
	$(CXX) $(CXX_FLAGS) -c edge.cxx

pseudograph.o: global.h pseudograph.h pseudograph.cxx
	$(CXX) $(CXX_FLAGS) -c pseudograph.cxx

logic_graph.o: global.h random.h edge.h schema.h pseudograph.h matrix.h binary_matrix.h graph.h logic_graph.h proposition.h propositional_system.h logic_graph.cxx
	$(CXX) $(CXX_FLAGS) -c logic_graph.cxx

polynomial_wrapper.o: global.h rational.h polynomial.h random.h polynomial.cxx polynomial_wrapper.cxx
	$(CXX) $(CXX_FLAGS) -c polynomial_wrapper.cxx

integer_polynomial_wrapper.o: global.h integer_polynomial.h integer_polynomial.cxx integer_polynomial_wrapper.cxx
	$(CXX) $(CXX_FLAGS) -c integer_polynomial_wrapper.cxx

variety_wrapper.o: global.h rational.h variety.h random.h variety.cxx variety_wrapper.cxx
	$(CXX) $(CXX_FLAGS) -c variety_wrapper.cxx

functional_equation_wrapper.o: global.h rational.h variety.h polynomial.h functional_equation.h functional_equation.cxx functional_equation_wrapper.cxx
	$(CXX) $(CXX_FLAGS) -c functional_equation_wrapper.cxx

propositional_system.o: global.h graph.h proposition.h propositional_system.h propositional_system.cxx binary_matrix.h directed_graph.h
	$(CXX) $(CXX_FLAGS) -c propositional_system.cxx

multitime.o: global.h multitime.h multitime.cxx
	$(CXX) $(CXX_FLAGS) -c multitime.cxx

vertex_wrapper.o: global.h vertex.h vertex.cxx vertex_wrapper.cxx
	$(CXX) $(CXX_FLAGS) -c vertex_wrapper.cxx

graph.o: graph.cxx graph.h edge.h schema.h pseudograph.h matrix.h binary_matrix.h random.h global.h
	$(CXX) $(CXX_FLAGS) -c graph.cxx

geometry_wrapper.o: geometry_wrapper.cxx geometry.cxx geometry.h matrix.h random.h global.h 
	$(CXX) $(CXX_FLAGS) -c geometry_wrapper.cxx

cell.o: cell.cxx cell.h global.h
	$(CXX) $(CXX_FLAGS) -c cell.cxx

nexus.o: nexus.cxx nexus.h cell.h schema.h graph.h edge.h random.h matrix.h binary_matrix.h pseudograph.h global.h
	$(CXX) $(CXX_FLAGS) -c nexus.cxx

group.o: group.cxx group.h word.h random.h global.h
	$(CXX) $(CXX_FLAGS) -c group.cxx

global.o: global.cxx random.h global.h
	$(CXX) $(CXX_FLAGS) -c global.cxx

random.o: random.cxx random.h global.h
	$(CXX) $(CXX_FLAGS) -c random.cxx

clean:
	rm -f $(OBJECTS)
	rm -f *~
	rm -f libsynarmosma.so

