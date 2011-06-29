/* Copyright (c) Zenna Tavares - zennatavares@gmail.com, 2010-2011 */
#pragma once

#include <algorithm>
#include <math.h>
#include <set>
#include <map>

namespace cliques {
/**
 @brief  Finds the variation of information between two partitions
 @param[in]  partition1     First partition
 @param[in]  partition2     Second Partition
 */
template<typename P>
double find_variation_of_information(P const &partition1, P const &partition2) {
	// Intersection of all sets
	std::map<int, std::set<int> > set_partition1, set_partition2;
	// number of nodes
	int num_elements = partition1.element_count();
	// map from comm_id to node_ids
	for (int i = 0; i < num_elements; ++i) {
		int comm_id = partition1.find_set(i);
		set_partition1[comm_id].insert(i);
	}
	for (int i = 0; i < num_elements; ++i) {
		int comm_id = partition2.find_set(i);
		set_partition2[comm_id].insert(i);
	}
	std::vector<int> inter_itr(num_elements);

	double var_inf = 0.0;
	for (auto set1 = set_partition1.begin(); set1
			!= set_partition1.end(); ++set1) {
		for (auto set2 = set_partition2.begin(); set2
				!= set_partition2.end(); ++set2) {

			auto end = std::set_intersection(set1->second.begin(),
					set1->second.end(), set2->second.begin(),
					set2->second.end(), inter_itr.begin());

			// number of elements both in p1:c_i and p2:c_j
			int n_ij = int(end - inter_itr.begin());

			if (n_ij != 0) {
				// number of elements in p1:ci
				int ni = set1->second.size();
				// number of elements in p2:cj
				int nj = set2->second.size();

				//std::cout << "ni " << ni << " nj " << nj << " nij " << n_ij << std::endl;

				double log_term = std::log( double(n_ij * n_ij)
						/ double(ni * nj));
				var_inf = var_inf + n_ij * log_term;
			}
		}
	}

	var_inf = -var_inf / (double(num_elements) * std::log(double(num_elements)));
	return var_inf;
}

}
