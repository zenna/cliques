// File: community.h
// -- community detection header file
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

#ifndef COMMUNITY_H
#define COMMUNITY_H

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <map>

#include "graph_binary.h"

#if defined(_WIN32)  || defined(WIN32)  || defined(__WIN32__) || defined(_WIN64)  || defined(WIN64) || defined(__WIN64__) || defined(_WINDOWS) || defined(__WINDOWS__)
#define __WIN__
#endif

#ifdef __WIN__ 
#include <windows.h>
#endif

using namespace std;

class Community {
 public:
     vector<double> neigh_weight;
  vector<unsigned int> neigh_pos;
  unsigned int neigh_last;

  // network to compute communities for
  Graph g;

  // nummber of nodes in the network and size of all vectors
  int size;

  // community to which each node belongs
  vector<int> n2c;

  // used to compute the modularity participation of each community
  vector<double> in;//,tot;

  double k_mean;

  // number of pass for one level computation
  // if -1, compute as many pass as needed to increase modularity
  int nb_pass;

  // Antoine
  // saves the initial number of nodes.
  int nb_nodes_init;



  // a new pass is computed if the last one has generated an increase 
  // greater than min_modularity
  // if 0. even a minor increase is enough to go for one more pass
  double min_modularity;
  
  double time;


  // constructors:
  // reads graph from file using graph constructor
  // type defined the weighted/unweighted status of the graph file
  Community (char *filename, int type, int nb_pass, double min_modularity, double t);
  // copy graph
  Community (Graph g, int nb_pass, double min_modularity, double t, int nb_nodes);
  
  Community (double * data, int length_data, int nbp, double minm, double timet, int type);
  
  ~Community();

  // display the community of each node
  void display();

  // Antoine
  // Modification of the function by adding int nb_nodes_per_comm_temp[]
  // containing the original number of nodes in each community.
  // remove the node from its current community with which it has dnodecomm links
  inline void remove(int node, int comm, double dnodecomm, int nb_nodes_per_comm_temp[]);

  // Antoine
  // Modification of the function by adding int nb_nodes_per_comm_temp[]
  // containing the original number of nodes in each community.
  // insert the node in comm with which it shares dnodecomm links
  inline void insert(int node, int comm, double dnodecomm, int nb_nodes_per_comm_temp[]);

  // compute the gain of modularity if node where inserted in comm
  // given that node has dnodecomm links to comm.  The formula is:
  // [(In(comm)+2d(node,comm))/2m - ((tot(comm)+deg(node))/2m)^2]-
  // [In(comm)/2m - (tot(comm)/2m)^2 - (deg(node)/2m)^2]
  // where In(comm)    = number of half-links strictly inside comm
  //       Tot(comm)   = number of half-links inside or outside comm (sum(degrees))
  //       d(node,com) = number of links from node to comm
  //       deg(node)   = node degree
  //       m           = number of links

  // Antoine
  // Modification of the function by adding int nb_nodes_per_comm_temp[]
  // containing the original number of nodes in each community.
  inline double modularity_gain(int node, int comm, double dnodecomm, int nb_nodes_per_comm_temp[]);

  // compute the set of neighboring communities of node
  // for each community, gives the number of links from node to comm
  void neigh_comm(unsigned int node);

  // compute the modularity of the current partition
  double modularity();

  // displays the graph of communities as computed by one_level
  void partition2graph();
  // displays the current partition (with communities renumbered from 0 to k-1)
  void display_partition(char *outfilename);
  
  
  vector<vector<int> >display_partition2(vector<vector<int> > output);
  
  // displays the current partition (with communities renumbered from 0 to k-1) on cerr
  void display_partitioncerr();

  // generates the binary graph of communities as computed by one_level
  Graph partition2graph_binary();

  // compute communities of the graph for one level
  // return the modularity
  bool one_level();
};

inline void
Community::remove(int node, int comm, double dnodecomm, int nb_nodes_per_comm_temp[]) {
  assert(node>=0 && node<size);

  // Antoine
  // You have to be careful not to forget suppressing
  // the nodes removed in the community where they were before
  // g.nb_nodes_per_comm[comm] -= nb_nodes_per_comm_temp[comm];
  // Only the nodes which were in this community BEFORE need to
  // be removed, not those that might have been added in the meantime
  // by the algorithm => necessary to have a temporary variable
  // nb_nodes_per_comm_temp.

   in[comm]  -= 2*dnodecomm + g.nb_selfloops(node);
  n2c[node]  = -1;
  g.nb_nodes_per_comm[comm] -= nb_nodes_per_comm_temp[node];
}

inline void
Community::insert(int node, int comm, double dnodecomm, int nb_nodes_per_comm_temp[]) {
  assert(node>=0 && node<size);

  // Antoine
  // You have to be careful not to forget to take into account
  // the nodes added in the community where they are now

  in[comm]  += 2*dnodecomm + g.nb_selfloops(node);
  n2c[node]=comm;
  g.nb_nodes_per_comm[comm] += nb_nodes_per_comm_temp[node];
}

// Antoine
// The gain in modularity is equal to (time/2m)*sum_ij(Aij) - #nodes_of_accepting_community*#nodes_of_giving_community*(1/N)^2
// Aij being all the links between the added nodes and the community where they are added
inline double
Community::modularity_gain(int node, int comm, double dnodecomm, int nb_nodes_per_comm_temp[]) {
  assert(node>=0 && node<size);
  double m2   = (double)g.total_weight;
  double dnc  = (double)dnodecomm;
  double gain = (time)*dnc-(g.nb_nodes_per_comm[comm]*nb_nodes_per_comm_temp[node])*(1.0/nb_nodes_init)*(1.0/nb_nodes_init)*m2;
  return gain;
}


#endif // COMMUNITY_H
