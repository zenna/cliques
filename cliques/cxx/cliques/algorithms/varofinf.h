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
	int num_elements = partition1.element_count();
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
	for (auto set1 = set_partition1.begin(); set1 != set_partition1.end(); ++set1) {
		for (auto set2 = set_partition2.begin(); set2 != set_partition2.end(); ++set2) {
			auto end = std::set_intersection(set1->second.begin(),
					set1->second.end(), set2->second.begin(),
					set2->second.end(), inter_itr.begin());
			int n_ij = int(end - inter_itr.begin());

			if (n_ij != 0) {
				int p1 = set1->second.size();
				int p2 = set2->second.size();
				double log_term = std::log(
						double(n_ij * n_ij) / double(p1 * p2));
				var_inf = var_inf + n_ij * log_term;
			}
		}
	}

	var_inf = -var_inf / double(num_elements);
	return var_inf;
}

}
