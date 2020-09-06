include Makefile.config

LD_FLAGS += -Wl,-rpath $(INSTALL_DIR)/lib

install:
	cd src; make
	mkdir -p $(INSTALL_DIR)/lib
	install -p src/libsynarmosma.so $(INSTALL_DIR)/lib/
	mkdir -p $(INSTALL_DIR)/include/synarmosma
	install -p -m 444 src/*.h $(INSTALL_DIR)/include/synarmosma/

test: install
	cd test; make

clean:
	cd src; make clean
	rm -f test/unit_test
	rm -f $(INSTALL_DIR)/lib/libsynarmosma.so
	rm -rf $(INSTALL_DIR)/include/synarmosma

