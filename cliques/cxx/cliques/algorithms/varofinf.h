/* Copyright (c) Zenna Tavares - zennatavares@gmail.com, 2010-2011 */
#ifndef CLIQUES_VAROFINF_H
#define CLIQUES_VAROFINF_H

#include <math.h>
#include <set>

namespace cliques
{
/**
@brief  Louvain method - greedy algorithm to find community structure of a network.
@param[in]  my_graph     graph to find partition of
@param[in]  quality_function     partition quality function object
*/
template <typename P>
float find_variation_of_information(P &partition1,P &partition2)
{
	float var_inf = 0.0;
    for (typename P::PartIterator p1itr = partition1.begin(); p1itr != partition1.end(); ++p1itr) {
    	for (typename P::PartIterator p2itr = partition2.begin(); p2itr != partition2.end(); ++p2itr) {
    		int p1 = partition1.set_size(p1itr);
    		int p2 = partition2.set_size(p2itr);
    		std::set<int> set_ij = P::intersection(p1itr,p2itr);
    		int n_ij = set_ij.size();

            if (n_ij != 0.0) {
                float log_term = std::log(float(n_ij*n_ij)/float(p1*p2));
                var_inf = var_inf + n_ij*log_term;
            }
    	}
    }
    var_inf = -var_inf/float(partition1.element_count());
    return var_inf;
}
}

#endif
