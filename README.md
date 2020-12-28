The Synarmosma library consists of approximately 21,000 lines of C++ code and has been developed
and tested in a Linux (x86, x86-64 and ARMv6) environment. The library's external dependencies include 
a modern C++ compiler as well as the [Boost library](https://www.boost.org), Victor Shoup's number-theoretic 
library [NTL](https://www.shoup.net/ntl), the GNU [multiprecision library](https://www.gmplib.org) and an 
implementation of LAPACK for certain dense matrix operations. The Synarmosma library is released under the 
GNU Public License version 3.0; see the <code>LICENSE.txt</code> file or the [Free Software Foundation](https://www.fsf.org/licensing) 
for more details.   

To compile the Synarmosma library, simply type the command <code>make</code> from the current directory, after 
which you can use the command <code>make test</code> to carry out some basic tests of certain classes and methods 
of the Synarmosma library. The command <code>make install</code> will install the library in the location which 
is specified by <code>INSTALL_DIR</code> in the <code>Makefile.config</code> file. This same file contains other 
parameters for controlling the compilation which you may modify if needed, for example to specify the location of 
the Boost or NTL library on your computer. The compilation of the Synarmosma library has been tested using both 
the GNU and Intel C++ compilers in a Linux environment. In principle it should also run under Windows using the 
Cygwin environment. Finally, you can use <code>make clean</code> to delete all the object files.   

Normally the only parameters that a user should need to modify in order to build the Synarmosma library
are located in the <code>Makefile.config</code> file, in which various compiler arguments and libraries are
specified. The build variables are all described in this file and in general the default values should be 
safe on most Unix-like platforms. Certain methods within the library have been parallelized using OpenMP,
namely for the <code>Geometry</code>, <code>Graph</code> and <code>Nexus</code> classes along with the 
library's matrix classes (<code>Matrix</code>, <code>Binary_Matrix</code>, <code>Integer_Matrix</code>), 
so that a C++ compiler able to support OpenMP is also desirable. The library's source code is documented 
using [Doxygen](https://www.doxygen.nl) with the configuration file <code>documentation/docs.config</code> 
that can be modified if desired. It is however sufficient to use the command <code>make docs</code> 
to create a set of HTML documents that will be installed in the <code>documentation/synarmosma</code> directory. 

For any questions, comments or suggestions, please contact <info@synarmosma.org>
