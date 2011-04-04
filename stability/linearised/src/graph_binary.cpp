// File: graph_binary.cpp
// -- graph handling source
//-----------------------------------------------------------------------------
// Community detection 
// Based on the article "Fast unfolding of community hierarchies in large networks"
// Copyright (C) 2008 V. Blondel, J.-L. Guillaume, R. Lambiotte, E. Lefebvre
//
// This program must not be distributed without agreement of the above mentionned authors.
//-----------------------------------------------------------------------------
// Author   : E. Lefebvre, adapted by J.-L. Guillaume
// Email    : jean-loup.guillaume@lip6.fr
// Location : Paris, France
// Time	    : February 2008
//-----------------------------------------------------------------------------
// see readme.txt for more details

#include <sys/mman.h>
#include "graph_binary.h"
#include <string.h>

Graph::Graph() : total_weight(0) {
  //nb_nodes     = 0;
  //nb_links     = 0;
  //total_weight = 0;
}

Graph::Graph(const char *filename, const bool weighted) {
	std::ifstream finput;
  finput.open(filename,std::fstream::in | std::fstream::binary);
  char curr_line[64];
	
	// check magic number
  finput.getline(curr_line,64,'\0');
	assert(!strcmp(curr_line,"MAG"));
	
	// check version number
	finput.getline(curr_line,64,'\0');
	assert(!strcmp(curr_line,"version"));
	finput.getline(curr_line,64,'\0');
	assert(!strcmp(curr_line,"1.0"));
	
	// check size of integer
	finput.getline(curr_line,64,'\0');
	assert(!strcmp(curr_line,"integer_size"));
	finput.getline(curr_line,64,'\0');
	assert(atoi(curr_line)==sizeof(int));

	// check size of double
	finput.getline(curr_line,64,'\0');
	assert(!strcmp(curr_line,"double_size"));
	finput.getline(curr_line,64,'\0');
	assert(atoi(curr_line)==sizeof(double));

  // read number of nodes on 4 bytes
	finput.getline(curr_line,64,'\0');
	assert(!strcmp(curr_line,"number_of_nodes"));
	int nb_nodes;
  finput.read((char *)&nb_nodes, sizeof(int));
  assert(nb_nodes>0);

  // read cumulative degree sequence: 4 bytes for each node
  // cum_degree[0]=degree(0); cum_degree[1]=degree(0)+degree(1), etc.
	finput.getline(curr_line,64,'\0');
	assert(!strcmp(curr_line,"start_nodes"));

	degrees.resize(nb_nodes);
  finput.read((char *)&degrees[0], nb_nodes*sizeof(int));

  // read links: one int per link (each link is counted twice)
  finput.getline(curr_line,64,'\0');
	assert(!strcmp(curr_line,"start_links"));
	int size_link_array=degrees[nb_nodes-1];
	links.resize(size_link_array);
  finput.read((char *)&links[0], size_link_array*sizeof(int));
  //cerr << "total : " << nb_links << endl;

  // IF weighted: read weights: 4 bytes for each link (each link is counted twice)
  if (weighted) {
		finput.getline(curr_line,64,'\0');
		assert(!strcmp(curr_line,"start_weights"));
		weights.resize(size_link_array);
    finput.read((char *)&weights[0], size_link_array*sizeof(int));
    total_weight=0;
    for (int i = 0 ; i<size_link_array ; i++) {
      total_weight += weights[i];
    }
  } else {
    //weights = NULL;
		// use instead for this test weights.size()==0
    total_weight = size_link_array;
  }
}

// generates a random graph using the benchmark approach
Graph::Graph(const int n1, const int k1, const int n2, const int k2, const int n3, const int k3) {
  srand(time(NULL));
  int nb_nodes = n1*n2*n3;

	std::vector<std::vector<int> > gr(nb_nodes);
  for (int i=0 ; i<nb_nodes ; i++)
    gr[i].resize(nb_nodes,0);

//  cerr << (k1*1.)/(n1*1.) << " " 
//       << (k2*1.)/(n1*n2*1.) << " " 
//       << (k3*1.)/(n1*n2*n3*1.) << endl;

  int nb_links = 0;
  for (int i=0 ; i<nb_nodes ; i++)
    for (int j=i+1 ; j<nb_nodes ; j++) {
      double v = rand()*1./RAND_MAX;
      if (i/n1==j/n1) { // i and j in the same subgroup
//	cout << i << " " << j << " 1 : " << v << " " << (k1*1.)/(n1-1.) ;
	if (v<=(k1*1.)/(n1-1.)) {
	  gr[i][j]=gr[j][i]=1;
	  nb_links++;
//	  cout << " : ok" ;
	}
//	cout << endl;
      } else if (i/(n1*n2)==j/(n1*n2)) { // i and j in the same group
//	cout << i << " " << j << " 2 : " << v << " " << (k2*1.)/(n1*(n2-1.)) ;
	if (v<=(k2*1.)/(n1*(n2-1.))) {
	  gr[i][j]=gr[j][i]=1;
	  nb_links++;
//	  cout << " : ok" ;
	}
//	cout << endl;
      } else { // i and j in different groups
//	cout << i << " " << j << " 3 : " << v << " " << (k3*1.)/(n1*n2*(n3-1.)) ;
	if (v<=(k3*1.)/(n1*n2*(n3-1.))) {
	  gr[i][j]=gr[j][i]=1;
	  nb_links++;
//	  cout << " : ok" ;
	}
//	cout << endl;
      }
    }

//  cerr << nb_links << endl;

  int size_link_array;
	size_link_array = total_weight = 2*nb_links;

	degrees.reserve(nb_nodes);
  for (int i=0 ; i<nb_nodes ; i++) {
    int d = 0;
    for (int j=0 ; j<nb_nodes ; j++)
      d+=gr[i][j];
    degrees.push_back(d);
  }
  for (int i=1 ; i<nb_nodes ; i++)
    degrees[i]+=degrees[i-1];

	links.reserve(size_link_array);
  for (int i=0 ; i<nb_nodes ; i++)
    for (int j=0 ; j<nb_nodes ; j++)
      if (gr[i][j]==1)
	links.push_back(j);

//  for (int i=0 ; i<nb_nodes ; i++)
//    cerr << degrees[i] << " " ;
//  cerr << endl;
}

// generates a random graph using the benchmark approach
Graph::Graph(const int n1, const int k1, const int n2, const int k2) {
  srand(getpid());

//  srand(time(NULL));
  int nb_nodes = n1*n2;

	std::vector<std::vector<int> > gr(nb_nodes);
  for (int i=0 ; i<nb_nodes ; i++)
    gr[i].resize(nb_nodes,0);

  int nb_links = 0;
  for (int i=0 ; i<nb_nodes ; i++)
    for (int j=i+1 ; j<nb_nodes ; j++) {
      double v = rand()*1./RAND_MAX;
      if (i/n1==j/n1) { // i and j in the same subgroup
	if (v<=(k1*1.)/(n1-1.)) {
	  gr[i][j]=gr[j][i]=1;
	  nb_links++;
	}
      } else { // i and j in different groups
	if (v<=(k2*1.)/(n1*(n2-1.))) {
	  gr[i][j]=gr[j][i]=1;
	  nb_links++;
	}
      }
    }
  int size_link_array;
	size_link_array = total_weight = 2*nb_links;

	degrees.reserve(nb_nodes);
  for (int i=0 ; i<nb_nodes ; i++) {
    int d = 0;
    for (int j=0 ; j<nb_nodes ; j++)
      d+=gr[i][j];
    degrees.push_back(d);
  }
  for (int i=1 ; i<nb_nodes ; i++)
    degrees[i]+=degrees[i-1];

	links.reserve(size_link_array);
  for (int i=0 ; i<nb_nodes ; i++)
    for (int j=0 ; j<nb_nodes ; j++)
      if (gr[i][j]==1)
	links.push_back(j);
}

// Think about replacing memcpy with copy at some point; "copy" code below nonfunctional!
Graph::Graph(const int nb_nodes, const int nb_links, const int t, const int *degrees_array, const int *links_array, const int *weights_array)
: total_weight(t) {
  int size_link_array = 2*nb_links;
	degrees.resize(nb_nodes);
	memcpy(&degrees[0],degrees_array,nb_nodes*sizeof(int));
	//copy (degrees_array[0],degrees_array[nb_nodes],&degrees[0]);
	links.resize(size_link_array);
	memcpy(&links[0],links_array,size_link_array*sizeof(int));
	//copy (links_array[0],links_array[size_link_array],&links[0]);
  weights.resize(size_link_array);
	memcpy(&weights[0],weights_array,size_link_array*sizeof(int));
	//copy (weights_array[0],weights_array[size_link_array],&weights[0]);
}
/*
Graph::Graph(const Graph& g)
: total_weight(g.total_weight),degrees(g.degrees),
	links(g.links),weights(g.weights){}
*/
Graph::~Graph() {
  total_weight = 0;
}

void
Graph::display() const {
  for (unsigned int node=0 ; node<degrees.size() ; node++) {
		std::pair<unsigned int *,int *> p = neighbors(node);
    for (int i=0 ; i<nb_neighbors(node) ; i++) {
      if (weights.size())
	std::cout << node << " " << *(p.first+i) << " " << *(p.second+i) << std::endl;
      else {
				std::cout << (node+1) << " " << (*(p.first+i)+1) << std::endl;
	//	cout << (node) << " " << (*(p.first+i)) << endl;
      }
    }
  }
}

void
Graph::display_binary(const char *outfile) {
	std::ofstream foutput;
  foutput.open(outfile ,std::fstream::out | std::fstream::binary);

  foutput << "MAG" << '\0' // Magic Number
	        << "version" << '\0' << "1.0" << '\0'
          << "integer_size" << '\0' << sizeof(int) << '\0'
					<< "double_size" << '\0' << sizeof(double) << '\0'
					<< "number_of_nodes" << '\0';
	int nb_nodes = degrees.size();
  foutput.write((char *)(&nb_nodes),sizeof(int));
	foutput << "start_nodes" << '\0';
  foutput.write((char *)(&degrees[0]),nb_nodes*sizeof(int));
	foutput << "start_links" << '\0';
  foutput.write((char *)(&links[0]),links.size()*sizeof(int));
  foutput << "start_weights" << '\0';
	foutput.write((char *)(&weights[0]),weights.size()*sizeof(int));
	foutput << "end_file" << '\0';

}


