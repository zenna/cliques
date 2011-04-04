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


#include <fstream>
#include "graph_binary.h"
#include "math.h"
extern "C" {
#include <math.h>
#include "mex.h"
#include "matrix.h"
}

#if defined(macintosh) || defined(MACOS) || defined(_MACOS)  || defined(__APPLE__)  || defined(_MAC)   || defined(MAC)  || defined(mac)  || defined(MACINTOSH)
#define __MAC__
#endif


#if defined(_WIN32)  || defined(WIN32)  || defined(__WIN32__) || defined(_WIN64)  || defined(WIN64) || defined(__WIN64__) || defined(_WINDOWS) || defined(__WINDOWS__)
#define __WIN__
#endif

#if defined(linux)  || defined(__linux__)  || defined(_linux) || defined(LINUX)  || defined(_LINUX) || defined(_UNIX) || defined(__UNIX__) || defined(__gnu_linux__) || defined(__unix__) || defined(UNIX) || defined(unix) || defined(sparc)
#define __lin__
#endif


#ifdef __WIN__
#endif

#ifdef __lin__
#include <sys/mman.h>
#endif

#ifdef __MAC__
#include <sys/mman.h>
#endif

Graph::Graph() {
  nb_nodes     = 0;
  nb_links     = 0;
  total_weight = 0;
}

Graph::Graph(char *filename, int type) {
  ifstream finput;
  finput.open(filename,fstream::in | fstream::binary);

  // read number of nodes on 4 bytes
  finput.read((char *)&nb_nodes, 4);

  // read cumulative degree sequence: 4 bytes for each node
  // cum_degree[0]=degree(0); cum_degree[1]=degree(0)+degree(1), etc.

  degrees = (unsigned long *)malloc((long)nb_nodes*sizeof(long));
  finput.read((char *)degrees, (long)nb_nodes*sizeof(int));

  // read links: 4 bytes for each link (each link is counted twice)
  nb_links=degrees[nb_nodes-1]/2;
  links = (unsigned int *)malloc((long)nb_links*2*sizeof(int));
  finput.read((char *)links, (long)nb_links*2*sizeof(int));  
  //cerr << "total : " << nb_links << endl;

  // IF WEIGHTED : read weights: 4 bytes for each link (each link is counted twice)
  if (type==WEIGHTED) {
    weights = (float *)malloc((long)nb_links*8);
    finput.read((char *)weights, (long)nb_links*8);  
    total_weight=0;
    for (unsigned int i = 0 ; i<nb_links*2 ; i++) {
      total_weight += weights[i];
    }
  } else {
    weights = NULL;
    total_weight = 2*nb_links;
  }

  // Antoine
  // New attribute being the backup of the number of nodes per community
  nb_nodes_per_comm = (int *)malloc((long)nb_nodes*4);
  for(unsigned int i = 0; i<nb_nodes; i++)
	  nb_nodes_per_comm[i] = 1;

}

Graph::Graph(double * data, int length_data, int type){


    nb_nodes = int(data[length_data-1])+1;
    
     degrees = (unsigned long *)malloc((long)nb_nodes*sizeof(unsigned long));//new unsigned long[nb_nodes];
     
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
     

   nb_links= degrees[nb_nodes-1]/2;  

    
   links =  (unsigned int *)malloc((long)nb_links*2*sizeof(unsigned int));//new unsigned int[nb_links];
   for(int i=0;i<length_data;i++)
    {
      links[i]=int(data[length_data+i]);
    }   

    
// IF WEIGHTED : read weights: 4 bytes for each link (each link is counted twice)
  weights = NULL;
  total_weight=length_data;
  if (type==WEIGHTED) {
  weights = (float *)malloc((long)nb_links*2*sizeof(float));//new float[nb_links];
  total_weight = 0;
    for(unsigned int i=0;i<nb_links*2;i++)
    {
      //printf("weights[i] = %d",data[2*length_data+/home/apd09/Documentsi]);
      weights[i]=(float)data[2*length_data+i];
      total_weight = total_weight+weights[i];
      //mexCallMATLAB(0,&result,0,&arguments,"pause");
    }
  }
  
  
  // Antoine
  // New attribute being the backup of the number of nodes per community
  nb_nodes_per_comm = (int *)malloc((long)nb_nodes*sizeof(int));
  for(unsigned int i = 0; i<nb_nodes; i++)
	  nb_nodes_per_comm[i] = 1;

 }

/*
// generates a random graph using the benchmark approach
Graph::Graph(int n1, int k1, int n2, int k2, int n3, int k3) {
  srand(time(NULL));
  nb_nodes = n1*n2*n3;

  vector<vector<int> > gr(nb_nodes);
  for (int i=0 ; i<nb_nodes ; i++)
    gr[i].resize(nb_nodes,0);

//  cerr << (k1*1.)/(n1*1.) << " " 
//       << (k2*1.)/(n1*n2*1.) << " " 
//       << (k3*1.)/(n1*n2*n3*1.) << endl;

  nb_links = 0;
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
  
  total_weight = 2*nb_links;
  weights      = NULL;

  degrees = (unsigned long *)malloc((long)nb_nodes*sizeof(long));
  for (int i=0 ; i<nb_nodes ; i++) {
    int d = 0;
    for (int j=0 ; j<nb_nodes ; j++)
      d+=gr[i][j];
    degrees[i]=d;
  }
  for (int i=1 ; i<nb_nodes ; i++)
    degrees[i]+=degrees[i-1];

  links = (unsigned int *)malloc((long)nb_links*sizeof(int));
  int pos=0;
  for (int i=0 ; i<nb_nodes ; i++)
    for (int j=0 ; j<nb_nodes ; j++)
      if (gr[i][j]==1)
	links[pos++]=j;

//  for (int i=0 ; i<nb_nodes ; i++)
//    cerr << degrees[i] << " " ;
//  cerr << endl;
}

// generates a random graph using the benchmark approach
Graph::Graph(int n1, int k1, int n2, int k2) {
  srand(getpid());

//  srand(time(NULL));
  nb_nodes = n1*n2;

  vector<vector<int> > gr(nb_nodes);
  for (int i=0 ; i<nb_nodes ; i++)
    gr[i].resize(nb_nodes,0);

  nb_links = 0;
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

  total_weight = 2*nb_links;
  weights      = NULL;

  degrees = (unsigned long*)malloc((long)nb_nodes*sizeof(long));
  for (int i=0 ; i<nb_nodes ; i++) {
    int d = 0;
    for (int j=0 ; j<nb_nodes ; j++)
      d+=gr[i][j];
    degrees[i]=d;
  }
  for (int i=1 ; i<nb_nodes ; i++)
    degrees[i]+=degrees[i-1];

  links = (unsigned int *)malloc((long)nb_links*sizeof(long));
  int pos=0;
  for (int i=0 ; i<nb_nodes ; i++)
    for (int j=0 ; j<nb_nodes ; j++)
      if (gr[i][j]==1)
	links[pos++]=j;
}
 */
/*
Graph::Graph(int n, int m, int t, int *d, int *l, int *w) {
  nb_nodes     = n;
  nb_links     = m;
  total_weight = t;
  degrees      = d;
  links        = l;
  weights      = w;
}
*/
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
