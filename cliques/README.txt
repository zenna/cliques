C++ 'cliques' library in cxx/cliques/

DEPENDENCIES

lemon graph library
boost or new g++  version with c++0x support
graphviz libray for drawing (apt-get install libgraphviz-dev)
doxygen to build documentation (apt-get install doxygen)
cmake for building (apt-get install cmake)

HOW TO BUILD

create a build directory e.g. >> mkdir ~/cliques/
>> cd ~/cliques/
>> cmake REPOSITORY_ROOT
>> make

DOCUMENTATION

>> cd REPOSITORY_ROOT/doc
>> doxygen
The doxumentation home page will now be doc/html/index.html

GUIDELINES FOR EXTENSION

*keep within the 'cliques' namespaces, do not use 'using' keyword
*keep modular by using templates and function objects (see existing code for reference)
*Document using doxygen format (see existing code)

make ;./scripts/create_communities -G ~/repos/graph-codes/cliques/data/graphs/random.edj -x barbell_n6; python ~/repos/graph-codes/cliques/cliques/scripts/create_json.py -x ~/cliques/barbell_n6 -o ~/repos/graph-codes/cliques/cliques.js/interactive/js/data/barbell_n6.jso
