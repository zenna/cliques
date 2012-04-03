#Cliques network and optimisation analysis library

**Cliques** is a set of libraries in C++, python, and Javascript for analysis of
1) Complex networks of varying sizes.
2) Landscapes induced by optimisation processes.

### Dependencies

lemon graph library
boost or new g++  version with c++0x support
doxygen to build documentation (apt-get install doxygen)
cmake for building (apt-get install cmake)

### How to build
Make a build directory (preferably outside the source code directory)

	mkdir ~/cliques/
	cd ~/cliques/
	cmake REPOSITORY_ROOT
	make

Note you may have to set flags to use a newer C++ compiler, and also to use gfortran, e.g

    cmake REPOSITORY_ROOT -D UseFortran=True -D CMAKE_CXX_COMPILER=g++-4.2

### Using JSON viewer

Note: Chrome has strict permissions for reading files out of the local file
system.  The easiest solution is to load the viewer using a local static web server.
Any static file web server will work; for example you can run Python's built-in server:

    cd REPOSITORY_ROOT
    python -m SimpleHTTPServer 8888

Then open your browser at <http://localhost:8888/>

### Documentation

    cd REPOSITORY_ROOT/doc
    doxygen

The documentation home page will now be doc/html/index.html

### Guidelines for extension

*keep within the 'cliques' namespaces, do not use 'using' keyword
*keep modular by using templates and function objects (see existing code for reference)
*Document using doxygen format (see existing code)