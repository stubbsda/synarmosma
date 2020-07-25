The Synarmosma library consists of approximately 20,000 lines of C++ code and has been developed
and tested in a Linux (x86, x86-64 and ARMv6) environment. The library's external dependencies include 
a modern C++ compiler as well as the [Boost library](https://www.boost.org), Victor Shoup's number-theoretic 
library [NTL](https://www.shoup.net/ntl), the GNU [multiprecision library](https://www.gmplib.org) and an 
implementation of LAPACK for certain dense matrix operations. The Synarmosma library is released under the 
GNU Public License version 3.0; see the <code>LICENSE.txt</code> file or the [Free Software Foundation](https://www.fsf.org/licensing) 
for more details.   

To install the Synarmosma library, enter the directory <code>synarmosma</code>, where the command
<code>make</code> will carry out the build process and <code>make test</code> to carry out some basic
tests of certain classes and methods of the Synarmosma library. The command <code>sudo make install</code>
will install the library in its default location, <code>/usr/local/synarmosma</code>. If you do not have root access
or do not wish to install the library in a system directory, you can modify the value of <code>INSTALL_DIR</code>
at the beginning of the Makefile to some other directory, such as <code>$HOME/synarmosma</code>. The
compilation of the Synarmosma library has been tested using both the GNU and Intel C++ compilers in a
Linux environment. In principle it should also run under Windows using the Cygwin environment. Finally
you can use <code>make clean</code> to delete all the object files.   

Normally the only parameters that a user should need to modify in order to build the Synarmosma library
are located at the beginning of the Makefile, where various compiler arguments and libraries are all
specified. The build variables are all described in the Makefile and in general the default values should
be safe on most Unix-like platforms. Certain methods within the library have been parallelized using
OpenMP, namely for the <code>Geometry</code>, <code>Graph</code> and <code>Nexus</code> classes along 
with the library's matrix classes (<code>Matrix</code>, <code>Binary_Matrix</code>, <code>Integer_Matrix</code>), 
so that a C++ compiler able to support OpenMP is also desirable. The library's source code is
documented using [Doxygen](https://www.doxygen.nl) with the configuration file <code>docs.config</code> that can be
modified as needed.

For any questions, comments or suggestions, please contact <info@synarmosma.org>
