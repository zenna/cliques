// File: lin_norm_lap-louvain.cpp
// -- community detection source file
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

#include "lin_norm_lap-louvain.h"

CommunityFinder::CommunityFinder(const Graph &gc, const int nb_pass, const double min_modularity)
  :g(gc), nb_pass(nb_pass),min_modularity(min_modularity)
{}


Community::Community(const Graph &gc, double time)
: size(gc.degrees.size()),g(gc), time(time)
{
  n2c.resize(size);
  in.resize(size);
  tot.resize(size);

  for (unsigned int i=0 ; i<size ; i++) {
    n2c[i] = i;
    in[i]  = g.nb_selfloops(i);
    tot[i] = g.weighted_degree(i);
  }
}

std::string
Community::display_string() const {
	std::ostringstream ost;
  ost << std::endl << "<" ;
  for (unsigned int i=0 ; i<size ; i++)
    ost << " " << i << "/" << n2c[i] << "/" << in[i] << "/" << tot[i] ;
  ost << ">" << std::endl;
	return ost.str();
}

double
Community::modularity() const {
  double q  = 1.0 - time;
  double m2 = (double)g.total_weight;

  for (unsigned int i=0 ; i<size ; i++) {
    if (tot[i]>0)
      q +=  time * (double)in[i]/m2 - ((double)tot[i]/m2)*((double)tot[i]/m2);
  }
  return q;
}

std::map<int,int>
Community::neigh_comm(unsigned int node) const {
  std::map<int,int> res;
  std::pair<unsigned int *,int *> p = g.neighbors(node);

  int deg = g.nb_neighbors(node);

  res.insert(std::make_pair(n2c[node],0));

  bool unweighted = (g.weights.size()==0);
  for (int i=0 ; i<deg ; i++) {
    unsigned int neigh        = *(p.first+i);
    int neigh_comm   = n2c[neigh];
    int neigh_weight = unweighted?1:*(p.second+i);

    if (neigh!=node) {
      std::map<int,int>::iterator it = res.find(neigh_comm);
      if (it!=res.end())
        it->second+=neigh_weight;
      else
        res.insert(std::make_pair(neigh_comm,neigh_weight));
    }
  }
  
  return res;
}


void
Community::partition2graph() const {
  std::vector<int> renumber(size, -1);
  for (unsigned int node=0 ; node<size ; node++) {
    renumber[n2c[node]]++;
  }

  int final=0;
  for (unsigned int i=0 ; i<size ; i++)
    if (renumber[i]!=-1)
      renumber[i]=final++;


  for (unsigned int i=0 ; i<size ; i++) {
    std::pair<unsigned int *,int *> p = g.neighbors(i);

    int deg = g.nb_neighbors(i);
    for (int j=0 ; j<deg ; j++) {
      int neigh = *(p.first+j);
      std::cout << renumber[n2c[i]] << " " << renumber[n2c[neigh]] << std::endl;
    }
  }
}

std::string
Community::string_partition() const {
  std::vector<int> renumber(size, -1);
  for (unsigned int node=0 ; node<size ; node++) {
    renumber[n2c[node]]++;
  }

  int final=0;
  for (unsigned int i=0 ; i<size ; i++)
    if (renumber[i]!=-1)
      renumber[i]=final++;

	std::ostringstream ost;
  for (unsigned int i=0 ; i<size ; i++)
    ost << i << " " << renumber[n2c[i]] << '\n';
	return ost.str();
}

Graph
Community::partition2graph_binary() const {
  std::vector<int> renumber(size, -1);
  for (unsigned int node=0 ; node<size ; node++) {
    renumber[n2c[node]]++;
  }

  int final=0;
  for (unsigned int i=0 ; i<size ; i++)
    if (renumber[i]!=-1)
      renumber[i]=final++;

  // Compute communities
  std::vector<std::vector<int> > comm_nodes(final);
  for (unsigned int node=0 ; node<size ; node++) {
    comm_nodes[renumber[n2c[node]]].push_back(node);
  }

  int comm_deg = comm_nodes.size();
  // unweigthed to weighted
	Graph g2;
  //g2.nb_nodes = comm_nodes.size();
	g2.degrees.resize(comm_deg);
  // First allocate space for the number of links in the parent graph. Will decrease later. 
  int memory_allocation = g.links.size();
	g2.links.reserve(memory_allocation);
  g2.weights.reserve(memory_allocation);
	
  bool unweighted = (g.weights.size()==0);
  int where = 0;
  for (int comm=0 ; comm<comm_deg ; comm++) {
    std::map<int,int> m;
    std::map<int,int>::iterator it;

    int comm_size = comm_nodes[comm].size();
    for (int node=0 ; node<comm_size ; node++) {
      std::pair<unsigned int *,int *> p = g.neighbors(comm_nodes[comm][node]);
      int deg = g.nb_neighbors(comm_nodes[comm][node]);
      for (int i=0 ; i<deg ; i++) {
        int neigh        = *(p.first+i);
        int neigh_comm   = renumber[n2c[neigh]];
        //int neigh_weight = (g.weights.size()==0)?1:*(p.second+i);
        int neigh_weight = unweighted?1:*(p.second+i);

        it = m.find(neigh_comm);
        if (it==m.end())
          m.insert(std::make_pair(neigh_comm, neigh_weight));
        else
          it->second+=neigh_weight;
      }
    }

    g2.degrees[comm]=(comm==0)?m.size():g2.degrees[comm-1]+m.size();

    for (it = m.begin() ; it!=m.end() ; it++) {
      g2.total_weight  += it->second;
      g2.links.push_back(it->first);
      g2.weights.push_back(it->second);
      where++;
    }
  }

  //cerr << where << " " << g2.nb_links << endl;
	g2.links.resize(where);
  g2.weights.resize(where);
 // cerr << g2.nb_links << endl;  
  return g2;
}

double
Community::one_level(const double min_modularity, const int nb_pass) {
  bool improvement = false;
  int nb_pass_done = 0;
  double new_mod   = modularity();
  double cur_mod   = new_mod;

  // repeat while 
  //   there is an improvement of modularity
  //   or there is an improvement of modularity greater than a given epsilon 
  //   or a predefined number of pass have been done

  std::vector<int> random_order(size);
  for (unsigned int i=0 ; i<size ; i++)
    random_order[i]=i;
  for (unsigned int i=0 ; i<size-1 ; i++) {
    int rand_pos = rand()%(size-i)+i;
    int tmp      = random_order[i];
    random_order[i] = random_order[rand_pos];
    random_order[rand_pos] = tmp;
  }

  do {
    cur_mod = new_mod;
    improvement = false;
    nb_pass_done++;

    // for each node: remove the node from its community and insert it in the best community
    for (unsigned int node_tmp=0 ; node_tmp<size ; node_tmp++) {
   //   int node = node_tmp;
      int node = random_order[node_tmp];
//      if (node%1000000==0) {fprintf(stderr,"%d ",node); fflush(stderr);}
      int node_comm     = n2c[node];

      // computation of all neighboring communities of current node
      std::map<int,int> ncomm   = neigh_comm(node);
      // remove node from its current community
      remove(node, node_comm, ncomm.find(node_comm)->second);

      // compute the nearest community for node
      // default choice for future insertion is the former community
      int best_comm        = node_comm;
      int best_nblinks     = 0;//ncomm.find(node_comm)->second;
      double best_increase = 0.;//modularity_gain(node, best_comm, best_nblinks);
      for (std::map<int,int>::iterator it=ncomm.begin() ; it!=ncomm.end() ; it++) {
        double increase = modularity_gain(node, it->first, it->second);
        if (increase>best_increase) {
          best_comm     = it->first;
          best_nblinks  = it->second;
          best_increase = increase;
        }
      }

      // insert node in the nearest community
      //      cerr << "insert " << node << " in " << best_comm << " " << best_increase << endl;
      insert(node, best_comm, best_nblinks);
     
      if (best_comm!=node_comm)
        improvement=true;
    }

    new_mod = modularity();
    //cerr << "pass number " << nb_pass_done << " of " << nb_pass 
    //     << " : " << new_mod << " " << cur_mod << endl;

  } while (improvement && new_mod-cur_mod>min_modularity && nb_pass_done!=nb_pass);
  return new_mod;
}

CommunityFinder::result
CommunityFinder::find_partition(const double time,const int display_level) const {
  std::ostringstream stdout_buffer;
  std::ostringstream stderr_buffer;

  //Community c2 = c;
  Community c2(g,time);
  time_t time_begin, time_end;
	std::time(&time_begin);
    
  const int numberinitial=c2.numberofcommunities();
  
  double qual = c2.modularity();
  double new_qual = c2.one_level(min_modularity, nb_pass);

  if (display_level==-1)
    stdout_buffer << c2.string_partition();

  Graph g2 = c2.partition2graph_binary();
  
  int level=0;

  int numberofcommunities=c2.numberofcommunities();

  while(new_qual-qual>min_modularity) {
    qual=new_qual;
    Community c2(g2, time);
	
  numberofcommunities=c2.numberofcommunities();
    
    new_qual = c2.one_level(min_modularity, nb_pass);
    
    if (display_level==-1)
      stdout_buffer << c2.string_partition();
    
    g2 = c2.partition2graph_binary();
    level++;
    
    if (level==display_level)
      g2.display();
    
  }
	std::time(&time_end);

  stderr_buffer << time << " " << numberinitial << " " << new_qual << " "<< numberofcommunities << " " << (time_end-time_begin) << std::endl;

  result findings(stdout_buffer.str(),stderr_buffer.str(),new_qual);
  return findings;
}

