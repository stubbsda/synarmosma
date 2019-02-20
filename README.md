The Synarmosma library consists of approximately 20,000 lines of C++ code and has been tested on 
the Linux (x86, x86-64 and ARMv6), OS X (x86-64) and OpenBSD (x86 and x86-64) platforms. The 
library's external dependencies consist of the Boost library (www.boost.org), Victor Shoup's 
number-theoretic library NTL (www.shoup.net/ntl) and an implementation of LAPACK for certain dense 
matrix operations. The Synarmosma library is released under the Gnu Public License version 3.0 
(see LICENSE.txt or www.fsf.org/licensing for details).   

To install the Synarmosma library, enter the directory "synarmosma", where the command "make" will carry 
out the build process and "make test" to carry out some basic tests of certain classes and methods 
of the Synarmosma library. The command "sudo make install" will install the library in its default 
location, /usr/local/synarmosma. If you don't have root access or don't wish to install the library 
in a system directory, you can modify the value of INSTALL_DIR at the beginning of the Makefile to 
some other directory, such as $HOME/synarmosma. The compilation of the Synarmosma library has been 
tested on several Unix platforms, including OS X/x86-64, OpenBSD (x86 and x86-64) and Linux (x86, 
x86-64 and ARMv6), using both the Gnu (g++) and Intel (icpc) compilers. In principle it should also 
run under Windows using the Cygwin environment. Finally you can use "make clean" to delete all the 
object files.   

Normally the only parameters that a user should need to modify in order to build the Synarmosma library 
are located at the beginning of the Makefile, where various compiler arguments and libraries are all 
specified. The build variables are all described in the Makefile and in general the default values should 
be safe on most Unix-like platforms. Certain methods within the library have been parallelized using 
OpenMP, namely for the Geometry, Graph, Nexus and Matrix classes, so that a C++ compiler able to support 
OpenMP is also desirable. A sample main.cpp file illustrating the use of the Variety template class 
is also included in this directory. 

For any questions, comments or suggestions, please contact info@synarmosma.org
