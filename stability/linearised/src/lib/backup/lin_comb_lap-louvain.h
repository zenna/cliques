// File: lin_norm_lap-louvain.h
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

#ifndef LIN_COMB_LAP_LOUVAIN_H
#define LIN_COMB_LAP_LOUVAIN_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <map>
#include <sstream>

#include "../abstract_community_finder.h"
#include "../graph_binary.h"

//Below is needed for compiling a shared object without name mangling
extern "C" {

class HGraph : public Graph {
 public:
  std::vector<int> nb_nodes_per_comm;
  HGraph(const Graph &gc);
  HGraph();
};

class Community {
 public:
  // community to which each node belongs
  std::vector<int> n2c;
  // used to compute the modularity participation of each community
  std::vector<int> in; 
  //std::vector<int> tot;
  double k_mean;

  // nummber of nodes in the network and size of all vectors
  const unsigned int size; // (nb_nodes_init)
  int nb_nodes_init; // I'm not sure what the difference is?

  // network to compute communities for
  HGraph g;
 
  const double time;
 public: 
  Community(const Graph &gc, double time, int nb_nodes);
  Community(const HGraph &gc, double time, int nb_nodes);

  // remove the node from its current community with which it has dnodecomm links
  inline void remove(unsigned int node, int comm, int dnodecomm, std::vector<int>& nb_nodes_per_comm_temp);

  // insert the node in comm with which it shares dnodecomm links
  inline void insert(unsigned int node, int comm, int dnodecomm, std::vector<int>& nb_nodes_per_comm_temp);

  // compute the modularity of the current partition
  double modularity() const;

  // compute the gain of modularity if node where inserted in comm
  // given that node has dnodecomm links to comm.  The formula is:
  // [(In(comm)+2d(node,comm))/2m - ((tot(comm)+deg(node))/2m)^2]-
  // [In(comm)/2m - (tot(comm)/2m)^2 - (deg(node)/2m)^2]
  // where In(comm)    = number of half-links strictly inside comm
  //       Tot(comm)   = number of half-links inside or outside comm (sum(degrees))
  //       d(node,com) = number of links from node to comm
  //       deg(node)   = node degree
  //       m           = number of links
  inline double modularity_gain(unsigned int node, int comm, int dnodecomm, std::vector<int>& nb_nodes_per_comm_temp) const;

  // compute the set of neighboring communities of node
  // for each community, gives the number of links from node to comm
  std::map<int,int> neigh_comm(unsigned int node) const;

  // returns a string containing the current partition (with communities renumbered from 0 to k-1
  std::string string_partition() const;

  // generates the binary graph of communities as computed by one_level
  HGraph partition2graph_binary() const;

  // displays the graph of communities as computed by one_level
  void partition2graph() const;

  // format as a string the community of each node
  std::string display_string() const;

  // returns the "size" of the partition; i.e. the number of communities
  unsigned int numberofcommunities() const;

  // compute communities of the graph for one level
  // return the modularity
  double one_level(const double min_modularity, const int nb_pass) ;

};

class CommunityFinder : public AbstractCommunityFinder {
 private:
  Community c;

  // number of pass for one level computation
  // if -1, compute as many pass as needed to increase modularity
  const int nb_pass;

  // a new pass is computed if the last one has generated an increase 
  // greater than min_modularity
  // if 0. even a minor increase is enough to go for one more pass
  const double min_modularity;
  
  const double time;

 public:
  // constructor:
  // Initialize by passing a graph reference
  CommunityFinder (const Graph &g, const int nb_pass,const double min_modularity,const double t);

  // Run one_level() repeatedly to give for an answer a single best partition.
  // return vector<string>(2) with the partition as the first entry and the summary
  // as the second.
  result find_partition(const int display_level) const;

  // make a public member function called quality to be more generalizable
  // quality() == modularity()
  //double quality() const;
  // rendered unnecessary by the better encapsulation of one_level in find_partition

  AbstractCommunityFinder *maker(const Graph &g, const int nb_pass,const double min_modularity,const double t);
};

inline void
Community::remove(unsigned int node, int comm, int dnodecomm, std::vector<int>& nb_nodes_per_comm_temp) {
  assert(node>=0 && node<size);

  in[comm]  -= 2*dnodecomm + g.nb_selfloops(node);
  n2c[node]  = -1;
  g.nb_nodes_per_comm[comm] -= nb_nodes_per_comm_temp[node];
}

inline void
Community::insert(unsigned int node, int comm, int dnodecomm, std::vector<int>& nb_nodes_per_comm_temp) {
  assert(node>=0 && node<size);

  in[comm]  += 2*dnodecomm + g.nb_selfloops(node);
  n2c[node]=comm;
  g.nb_nodes_per_comm[comm] += nb_nodes_per_comm_temp[node];
}

inline double
Community::modularity_gain(unsigned int node, int comm, int dnodecomm, std::vector<int>& nb_nodes_per_comm_temp) const {
  assert(node>=0 && node<size);

  double m2   = (double)g.total_weight;
  double dnc  = (double)dnodecomm;
  double gain = (time)*dnc-(g.nb_nodes_per_comm[comm]*nb_nodes_per_comm_temp[node])*(1.0/nb_nodes_init)*(1.0/nb_nodes_init)*m2;
  return gain;
}

inline unsigned int
Community::numberofcommunities() const{return g.degrees.size();}

AbstractCommunityFinder *maker(const Graph &g, const int nb_pass,const double min_modularity,const double t) 
{
  //return new Community(const Graph &g, const int nb_pass,const double min_modularity,const double t);
  return new CommunityFinder(g, nb_pass, min_modularity, t);
}

class proxy {
  public:
    // register the maker with the factory
    //proxy(const Graph &g, const int nb_pass,const double min_modularity,const double t){
    proxy(){
      //factory["community"] = maker(g, nb_pass, min_modularity, t);
      factory["lin_comb_lap-louvain"] = maker;
    }
};

// Our one instance of the proxy
proxy p;

}


#endif // LIN_COMB_LAP_LOUVAIN_H
