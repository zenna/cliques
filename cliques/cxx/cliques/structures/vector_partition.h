/* Copyright (c) Z Tavares, M Schaub - 2010-2011 */
#ifndef VECTOR_PARTITION_H
#define VECTOR_PARTITION_H

#include <vector>
#include <map>

namespace cliques {
/**
 @brief  A partition where the position in the vector denotes the node_id
 and its value determines its set id.
 Good when you need constant time assignment of a node to a set
 and removal of a node from a set
 */
class VectorPartition {
private:
	int num_nodes;
	std::vector<int> partition_vector;

public:
	//#################### CONSTRUCTORS ####################

	// construct empty partition
	explicit VectorPartition(int num_nodes) :
		num_nodes(num_nodes), partition_vector(std::vector<int>(num_nodes, -1)) {
	}

	// construct partition from vector
	explicit VectorPartition(std::vector<int> partition) :
		partition_vector(partition) {
	}
	// what is that for a horribly over complicated syntax?? --M
	void initialise_as_singletons() {
		for (std::vector<int>::iterator itr = partition_vector.begin(); itr
				!= partition_vector.end(); ++itr) {
			*itr = itr - partition_vector.begin();

			double test = *itr;
		}
	}

	//#################### PUBLIC METHODS ####################
public:
	int find_set(int node_id) const {
		return partition_vector[node_id];
	}

	void isolate_node(int node_id) {
		partition_vector[node_id] = -1;
	}

	void add_node_to_set(int node_id, int set_id) {
		partition_vector[node_id] = set_id;
	}

	void normalise_ids() {
		int start_num = 0;
		std::map<int, int> set_to_node;
		std::map<int, int> set_old_to_new;
		for (std::vector<int>::iterator itr = partition_vector.begin(); itr
				!= partition_vector.end(); ++itr) {

			std::map<int, int>::iterator old_set = set_old_to_new.find(*itr);
			if (old_set == set_to_node.end()) {
				set_old_to_new[*itr] = start_num;
				*itr = start_num;
				start_num++;
			} else {
				*itr = old_set->second;
			}
		}
	}

	int element_count() {
		return partition_vector.size();
	}

	int set_count() {
		std::vector<int> count_vector(num_nodes, 0);
		for (std::vector<int>::iterator itr = partition_vector.begin(); itr
				!= partition_vector.end(); ++itr) {
			if (*itr == -1)
				// TODO warning: there are unassigned nodes
				return -1;
			else
				count_vector[*itr] = 1;
		}
		int count = 0;
		for (std::vector<int>::iterator itr = count_vector.begin(); itr
				!= count_vector.end(); ++itr) {
			count += *itr;
		}
		return count;
	}

	//#################### ITERATORS ####################
	class NodeIterator;

	class PartIterator {
		int alpha;
	};

	PartIterator begin() {
		PartIterator alpha;
		return alpha;
	}

	PartIterator end() {
		PartIterator alpha;
		return alpha;
	}

	class NodeIterator {
		int alpha;
	};

};

}

#endif
