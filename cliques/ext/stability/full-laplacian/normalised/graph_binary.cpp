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
#include <fstream>
#include "graph_binary.h"
#include "math.h"
extern "C" {
#include <math.h>
#include "mex.h"
}


Graph::Graph() {
  nb_nodes     = 0;
  nb_links     = 0;
  total_weight = 0;
}

Graph::Graph(char *filename, char *filename_w, int type) {
  ifstream finput;
  finput.open(filename,fstream::in | fstream::binary);

  // Read number of nodes on 4 bytes
  finput.read((char *)&nb_nodes, 4);
  assert(finput.rdstate() == ios::goodbit);

  // Read cumulative degree sequence: 8 bytes for each node
  // cum_degree[0]=degree(0); cum_degree[1]=degree(0)+degree(1), etc.

  // ligne modifiee par Antoine : 4 bytes au lieu de 8
  degrees = (unsigned long *)malloc((long)nb_nodes*4);
  // ligne modifiee par Antoine : 4 bytes au lieu de 8
  finput.read((char *)degrees, (long)nb_nodes*4);
  assert(finput.rdstate() == ios::goodbit);

  // Read links: 4 bytes for each link (each link is counted twice)
  nb_links=degrees[nb_nodes-1];

  // Ligne modifiee par Antoine : il faut mettre l'affichage apres la redefinition de nb_links
  //cerr << "total : " << nb_nodes << " " << nb_links << endl;


  links = (unsigned int *)malloc((long)nb_links*4);
  finput.read((char *)links, (long)nb_links*4);  
//  assert(finput.rdstate() == ios::goodbit);
  assert(links);

  // IF WEIGHTED : read weights: 4 bytes for each link (each link is counted twice)
  weights = NULL;
  total_weight=0;
  if (type==WEIGHTED) {
    ifstream finput_w;
    finput_w.open(filename_w,fstream::in | fstream::binary);

    weights = (float *)malloc((long)nb_links*4);
    finput_w.read((char *)weights, (long)nb_links*4);  
//    assert(finput_w.rdstate() == ios::goodbit);

//    assert(check_symmetry());
  }    

  // Compute total weight
  for (unsigned int i=0 ; i<nb_nodes ; i++) {
    total_weight += (double)weighted_degree(i);
  }

  //cerr << "total : " << nb_nodes << " " << nb_links << endl;

}

 Graph::Graph(double * data, int length_data){
         //Print the integer avg of each col to matlab console
     
    // mxArray *result;
   // mxArray *arguments;
   // printf("In Graph constructor.");
    
   // mexCallMATLAB(0,&result,0,&arguments,"pause");

    nb_nodes = int(data[length_data-1])+1;
  //  printf("nb_nodes = %d", nb_nodes);
    
    //mexCallMATLAB(0,&result,0,&arguments,"pause");
    
     degrees = (unsigned long *)malloc((long)nb_nodes*4);//new unsigned long[nb_nodes];
     
     int tmp_node = -1;
     int j = 0;
     int tot = 0;
   for(int i=0;i<length_data;i++)
    {
        tot += 1;
       if(int(data[i]) == tmp_node)
           degrees[tmp_node] += 1;
       else{
           tmp_node ++;
           degrees[tmp_node] = tot;
       }
   }
     
    // printf("degrees[0] = %d", degrees[0]);
    
   // mexCallMATLAB(0,&result,0,&arguments,"pause");
     
   nb_links= degrees[nb_nodes-1];  
   links =  (unsigned int *)malloc((long)nb_links*4);//new unsigned int[nb_links];
   for(int i=0;i<length_data;i++)
    {
      links[i]=int(data[length_data+i]);
    }   
   
   //printf("links[1] = %d", links[1]);
    
    //mexCallMATLAB(0,&result,0,&arguments,"pause");
    
  weights = (float *)malloc((long)nb_links*4);//new float[nb_links];
  total_weight = 0;
    for(int i=0;i<nb_links;i++)
    {
      weights[i]=(float)data[2*length_data+i];
      total_weight += (double)weights[i];
    }
  
    //printf("\n\nfunction data[1,2] = %f, data[1,3] = %f", links[0], weights[0]);
    
   // mexCallMATLAB(0,&result,0,&arguments,"pause");

 }

Graph::Graph(int n, int m, double t, int *d, int *l, float *w) {
/*  nb_nodes     = n;
  nb_links     = m;
  total_weight = t;
  degrees      = d;
  links        = l;
  weights      = w;*/
}


void
Graph::display() {
  for (unsigned int node=0 ; node<nb_nodes ; node++) {
    pair<unsigned int *, float *> p = neighbors(node);
    for (unsigned int i=0 ; i<nb_neighbors(node) ; i++) {
      if (node<=*(p.first+i)) {
	if (weights!=NULL)
	  cout << node << " " << *(p.first+i) << " " << *(p.second+i) << endl;
	else
	  cout << node << " " << *(p.first+i) << endl;
      }
    }   
  }
}

void
Graph::display_reverse() {
  for (unsigned int node=0 ; node<nb_nodes ; node++) {
    pair<unsigned int *, float *> p = neighbors(node);
    for (unsigned int i=0 ; i<nb_neighbors(node) ; i++) {
      if (node>*(p.first+i)) {
	if (weights!=NULL)
	  cout << *(p.first+i) << " " << node << " " << *(p.second+i) << endl;
	else
	  cout << *(p.first+i) << " " << node << endl;
      }
    }   
  }
}


bool
Graph::check_symmetry() {
  int error=0;
  for (unsigned int node=0 ; node<nb_nodes ; node++) {
    pair<unsigned int *, float *> p = neighbors(node);
    for (unsigned int i=0 ; i<nb_neighbors(node) ; i++) {
      unsigned int neigh = *(p.first+i);
      float weight = *(p.second+i);
      
      pair<unsigned int *, float *> p_neigh = neighbors(neigh);
      for (unsigned int j=0 ; j<nb_neighbors(neigh) ; j++) {
	unsigned int neigh_neigh = *(p_neigh.first+j);
	float neigh_weight = *(p_neigh.second+j);

	if (node==neigh_neigh && weight!=neigh_weight) {
	  cout << node << " " << neigh << " " << weight << " " << neigh_weight << endl;
	  if (error++==10)
	    exit(0);
	}
      }
    }
  }
  return (error==0);
}


void
Graph::display_binary(char *outfile) {
  ofstream foutput;
  foutput.open(outfile ,fstream::out | fstream::binary);

  foutput.write((char *)(&nb_nodes),4);
  foutput.write((char *)(degrees),4*nb_nodes);
  foutput.write((char *)(links),8*nb_links);
}
