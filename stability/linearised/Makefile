#!/bin/bash

CXX=g++
CXXFLAGS= -ansi -O3 -Wall -fPIC 
LDFLAGS= -ansi -lm -Wall -fPIC
LIBRARIES=$(patsubst src/lib/%.cpp,lib/lib%.so,$(wildcard src/lib/*.cpp))
EXEC=community convert hierarchy random_graph bin2txt $(LIBRARIES)
OBJ1= graph_binary.o community.o
OBJ2= graph.o
LIBS = -ldl
VPATH = src

all: $(EXEC)

#community : $(OBJ1) main_community.o -lboost_program_options
community : graph_binary.o main_community.o 
	$(CXX) -o $@ $^ $(LDFLAGS) -rdynamic -lboost_program_options $(LIBS)

lib/lib%.so : lib/%.o
	$(CXX) -shared -Wl,-soname,lib$*.so -o lib/lib$*.so lib/$*.o $(LDFLAGS)
 
convert : $(OBJ2) main_convert.o
	$(CXX) -o $@ $^ $(LDFLAGS)

hierarchy : main_hierarchy.o
	$(CXX) -o $@ $^ $(LDFLAGS)

random_graph : graph_binary.o main_random_graph.o
	$(CXX) -o $@ $^ $(LDFLAGS) -lboost_program_options 

bin2txt : graph_binary.o main_bin2txt.o
	$(CXX) -o $@ $^ $(LDFLAGS) -lboost_program_options 

##########################################
# Generic rules
##########################################

%.o: %.cpp %.h
	$(CXX) -o $@ -c $< $(CXXFLAGS)

%.o: %.cpp
	$(CXX) -o $@ -c $< $(CXXFLAGS)

lib/%.o : lib/%.cpp lib/%.h
	$(CXX) -o $@ -c $< $(CXXFLAGS) 

.SECONDARY: $(wildcard lib/*.o)

clean:
	rm -f *.o *~ lib/*.o lib/*~ $(EXEC)
