// File: abstract_community_finder.h
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

#ifndef ABSTRACT_COMMUNITY_FINDER_H
#define ABSTRACT_COMMUNITY_FINDER_H

//#include <stdlib>
//#include <stdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <map>
#include <sstream>
#include <string>

#include "graph_binary.h"


extern "C" {

//using namespace std;

class AbstractCommunityFinder {
 public:
  virtual ~AbstractCommunityFinder () {};

  // structure for the result of findpartition()
  struct result{
    std::string partition;
    std::string partition_summary;
    double quality;

    result(std::string part, std::string part_sum, double qual)
    :partition(part),partition_summary(part_sum),quality(qual) {}
  };

  virtual result find_partition(const double time, const int display_level) const = 0;
};


// typedef to make it easier to set up our factory
typedef AbstractCommunityFinder *maker_t(const Graph &g, const int nb_pass,const double min_modularity);
// our global factory
extern std::map<std::string, maker_t *, std::less<std::string> > factory;

}

#endif // ABSTRACT_COMMUNITY_H
