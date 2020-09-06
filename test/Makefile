include ../Makefile.config

LD_FLAGS += -Wl,-rpath $(INSTALL_DIR)/lib

test:
	$(CXX) -std=c++14 -Wall -march=native -Wl,-rpath $(INSTALL_DIR)/lib -I$(INSTALL_DIR)/include -L$(INSTALL_DIR)/lib -o unit_test unit_testing.cxx -lsynarmosma $(LIBS)
	./test.sh
	rm unit_test
