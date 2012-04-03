/* Copyright (c) Z Tavares, M Schaub - 2010-2011 */
#ifndef VECTOR_PARTITION_H
#define VECTOR_PARTITION_H

#include <vector>
#include <set>
#include <map>
#include <iostream>

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
	bool is_normalised;

public:
	//#################### CONSTRUCTORS ####################

	// construct empty partition
	explicit VectorPartition(int num_nodes) :
		num_nodes(num_nodes),
				partition_vector(std::vector<int>(num_nodes, -1)),
				is_normalised(false) {
	}

	// construct partition with initial_set
	explicit VectorPartition(int num_nodes, int initial_set) :
		num_nodes(num_nodes), partition_vector(std::vector<int>(num_nodes,
				initial_set)), is_normalised(false) {
	}

	// construct partition from vector
	explicit VectorPartition(std::vector<int> partition) :
		num_nodes(partition.size()), partition_vector(partition),
				is_normalised(false) {
	}
	// what is that for a horribly over complicated syntax?? --M
	void initialise_as_singletons() {
		for (std::vector<int>::iterator itr = partition_vector.begin(); itr
				!= partition_vector.end(); ++itr) {
			*itr = itr - partition_vector.begin();
		}
		is_normalised = true;
	}

	void initialise_as_global() {
		partition_vector = std::vector<int>(partition_vector.size(), 0);
	}

	//#################### PUBLIC METHODS ####################
public:
	int find_set(int node_id) const {
		return partition_vector[node_id];
	}

	void isolate_node(int node_id) {
		partition_vector[node_id] = -1;
		is_normalised = false;
	}

	void add_node_to_set(int node_id, int set_id) {
		partition_vector[node_id] = set_id;
		is_normalised = false;
	}

	std::map<int, int> normalise_ids() {
		std::map<int, int> set_new_to_old;
		if (!is_normalised) {
			int start_num = 0;
			std::map<int, int> set_old_to_new;
			for (std::vector<int>::iterator itr = partition_vector.begin(); itr
					!= partition_vector.end(); ++itr) {

				std::map<int, int>::iterator old_set =
						set_old_to_new.find(*itr);
				if (old_set == set_old_to_new.end()) {
					set_old_to_new[*itr] = start_num;
					set_new_to_old[start_num] = *itr;
					*itr = start_num;
					start_num++;
				} else {
					*itr = old_set->second;
				}
			}
			is_normalised = true;
			return set_new_to_old;
		} else {
			for (std::vector<int>::iterator itr = partition_vector.begin(); itr
					!= partition_vector.end(); ++itr) {
				set_new_to_old[*itr] = *itr;
			}
			return set_new_to_old;
		}
	}

	int element_count() const {
		return partition_vector.size();
	}

	int set_count() {
		std::set<int> seen_nodes;
		for (std::vector<int>::iterator itr = partition_vector.begin(); itr
				!= partition_vector.end(); ++itr) {
			if (*itr != -1) {
				seen_nodes.insert(*itr);
			}
		}
		return seen_nodes.size();
	}

	std::vector<int> get_nodes_from_set(int set_id) {
		std::vector<int> nodes_in_set;
		for (int i = 0; i < num_nodes; ++i) {
			if (partition_vector[i] == set_id) {
				nodes_in_set.push_back(i);
			}
		}
		return nodes_in_set;
	}

	std::vector<int> return_partition_vector() {
		return partition_vector;
	}

	//#################### Operators ####################
	bool operator==(const VectorPartition& other) const {
		if (this->is_normalised && other.is_normalised) {
			return (other.partition_vector == this->partition_vector);
		} else {
			VectorPartition a(*this);
			VectorPartition b(other);

			a.normalise_ids();
			b.normalise_ids();
			return a.partition_vector == b.partition_vector;
		}
		/*if (partition_vector.size() != other.partition_vector.size()) {
		 return false;
		 }
		 const int MAGIC_NUM = -321;
		 int num_elements = partition_vector.size();

		 std::vector<int> this_to_other(num_elements, MAGIC_NUM);
		 std::vector<int> other_to_this(num_elements, MAGIC_NUM);

		 for (int i=0; i<num_elements; ++i) {
		 int a = partition_vector[i];
		 int b = other.partition_vector[i];

		 // TODO:Replace with exception handling or generalise to any set num
		 if (a >= num_elements || b >= num_elements || a < 0 || b < 0) {
		 std::cout << "out of bounds!";
		 }

		 if (this_to_other[a] == MAGIC_NUM && other_to_this[b] == MAGIC_NUM) {
		 this_to_other[a] = b;
		 other_to_this[b] = a;
		 }
		 else if (this_to_other[a] != b || other_to_this[b] != a) {
		 return false;
		 }
		 }
		 return true;*/
	}

	bool operator!=(const VectorPartition& other) const {
		return !(partition_vector == other.partition_vector);
	}

	//#################### ITERATOR SETUP ####################
	class PartIterator;
	class NodeIterator;

	PartIterator begin() {
		std::vector<int>::iterator tmp_part_itr = partition_vector.begin();
		return PartIterator(tmp_part_itr, *this);
	}

	PartIterator end() {
		std::vector<int>::iterator tmp_part_itr = partition_vector.end();
		return PartIterator(tmp_part_itr, *this);
	}

	/*NodeIterator start() {
	 std::vector<int>::iterator tmp_part_itr =
	 partition_vector.begin();
	 return NodeIterator(tmp_part_itr, *this);
	 }

	 NodeIterator finish() {
	 std::vector<int>::iterator tmp_part_itr =
	 partition_vector.end();
	 return NodeIterator(tmp_part_itr, *this);
	 }*/

	//#################### ITERATORS ####################
	class PartIterator {
	private:
		std::set<int> visited_sets;
		std::vector<int>::iterator part_itr;
		VectorPartition &parent_partition;

	public:
		PartIterator(std::vector<int>::iterator part_itr,
				VectorPartition &parent_partition) :
			part_itr(part_itr), parent_partition(parent_partition) {
		}

		PartIterator& operator=(const PartIterator& other) {
			part_itr = other.part_itr;
			visited_sets = other.visited_sets;
			return (*this);
		}

		bool operator==(const PartIterator& other) {
			return (part_itr == other.part_itr);
		}

		bool operator!=(const PartIterator& other) {
			return (part_itr != other.part_itr);
		}

		int operator*() {
			return (*part_itr);
		}

		PartIterator& operator++() {
			int current_set = *part_itr;
			visited_sets.insert(current_set);

			if (visited_sets.size() == parent_partition.partition_vector.size()) {
				part_itr = parent_partition.partition_vector.end();
			} else {
				do {
					++part_itr;
				} while (visited_sets.find(*part_itr) != visited_sets.end());
			}
			return *this;
		}
	};

	/*class NodeIterator {
	 NodeIterator(std::vector<int>::iterator part_itr,
	 VectorPartition &parent_partition) :
	 part_itr(part_itr), parent_partition(parent_partition) {}

	 NodeIterator& operator=(const PartIterator& other) {
	 part_itr = other.part_itr;
	 visited_sets = other.visited_sets;
	 return (*this);
	 }

	 bool operator==(const NodeIterator& other) {
	 return (part_itr == other.part_itr);
	 }

	 bool operator!=(const NodeIterator& other) {
	 return (part_itr != other.part_itr);
	 }

	 int operator*() {
	 return (*part_itr);
	 }

	 PartIterator& operator++() {
	 current_node++;
	 return *this;
	 }
	 };*/

};

}

#endif
