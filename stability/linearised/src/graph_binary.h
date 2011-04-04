// File: graph_binary.h
// -- graph handling header file
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

#ifndef GRAPH_H
#define GRAPH_H

//#include <stdlib.h>
//#include <stdio.h>
#include <assert.h>
//#include <malloc.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <map>
#include <algorithm>

//using namespace std;

extern "C" {

class Graph {
//  bool weighted;
 public:
	//nb_nodes = degrees.size()
  int total_weight;  

	std::vector<int> degrees;
	std::vector<int> links;
	std::vector<int> weights;

  Graph();

  // binary file format is
  // 4 bytes for the number of nodes in the graph
  // 4*(nb_nodes) bytes for the cumulative degree for each node:
  //    deg(0)=degrees[0]
  //    deg(k)=degrees[k]-degrees[k-1]
  // 4*(sum_degrees) bytes for the links
  // IF weighted 4*(sum_degrees) bytes for the weights
  //Graph(char *filename, bool weighted);
  Graph(const char *filename,const bool weighted);

  Graph(const int nb_nodes,const int nb_links,const int total_weight,const int *degrees,const int *links,const int *weights);
  
  Graph(const int n1,const int k1,const int n2,const int k2,const int n3,const int k3);
  Graph(const int n1,const int k1,const int n2,const int k2);
//	Graph(const Graph& g);

  ~Graph();
  
  void display(void) const;
  void display_binary(const char *outfile);

  // return the number of neighbors (degree) of the node
  inline int nb_neighbors(const unsigned int node) const;

  // return the number of self loops of the node
  inline int nb_selfloops(const unsigned int node) const;

  // return the weighted degree of the node
  inline int weighted_degree(const unsigned int node) const;

  // return pointers to the first neighbor and first weight of the node
  inline std::pair<unsigned int *,int *> neighbors(unsigned int node) const;
};


inline int
Graph::nb_neighbors(unsigned int node) const {
  assert(node>=0 && node<degrees.size());
  if (node==0)
    return degrees[0];
  else
    return degrees[node]-degrees[node-1];
}

inline int
Graph::nb_selfloops(unsigned int node) const {
  assert(node>=0 && node<degrees.size());

	std::pair< const unsigned int *, const int *> p = neighbors(node);
  for (int i=0 ; i<nb_neighbors(node) ; i++) {
    if (*(p.first+i)==node) {
      if (weights.size())
	return *(p.second+i);
      else 
	return 1;
    }
  }
  return 0;
}

inline int
Graph::weighted_degree(unsigned int node) const {
   assert(node>=0 && node<degrees.size());

	 std::pair<unsigned int *,int *> p = neighbors(node);
   if (p.second==NULL)
     return nb_neighbors(node);
   else {
     int res = 0;
     for (int i=0 ; i<nb_neighbors(node) ; i++)
       res += *(p.second+i);
     return res;
   }
}

inline std::pair< unsigned int *, int *>
Graph::neighbors(unsigned int node) const {
	//cerr << node << " < " << degrees.size() << endl;
  assert(node>=0 && node<degrees.size());

	//const int * links_pointer = &links[0];
	//const int * weights_pointer = &weights[0];
  if (node==0)
    return std::make_pair( (unsigned int *)&links[0], (int *)&weights[0]);
  else if (weights.size())
    return std::make_pair( (unsigned int *)&links[degrees[node-1]],  (int *)&weights[degrees[node-1]]);
  else
    return std::make_pair( (unsigned int *)&links[degrees[node-1]], (int *)0);
}


}
#endif // GRAPH_H
