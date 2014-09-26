/* Copyright (c) Z Tavares - 2010-2011 */
#ifndef TINY_PARTITION_H
#define TINY_PARTITION_H

#include <cmath>
#include <vector>
#include <map>
#include <iostream>

namespace clq {
/**
 @brief  A partition where the position in the vector denotes the node_id
 and its value determines its set id.

 Good when you need space efficiency and constant (albeit probably slower
 than vector) time assignment of a node to a set and removal of a node from a
 set
 */
class TinyPartition {
public:
	int num_nodes;
	unsigned char num_chars;
	unsigned char* partition_array;

public:
	//#################### CONSTRUCTORS ####################

	// construct empty partition
	explicit TinyPartition(int num_nodes) :
	num_nodes(num_nodes) {
		bits_per_node = ceil(float(std::log(num_nodes) / std::log(2)));
		float total_bits_required = bits_per_node * num_node;
		int num_chars_required = total_bits_required / (sizeof(unsigned char) * 8);
		partition_array = new unsigned char(num_chars_required);
	}

	~TinyPartition() {
		delete[] partition_array;
	}

	// what is that for a horribly over complicated syntax?? --M
	void initialise_as_singletons() {
		for (int i=0; i<num_nodes; ++i) {
			add_node_to_set(i,i);
		}
	}

private:
	void left_shift() {
		unsigned char length = 10;
		unsigned char data[10] = {0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0A,0xBC};
		unsigned char *shift = data;
		while (shift < data+(length-2)) {
			*shift = (*(shift+1)&0x0F)<<4 | (*(shift+2)&0xF0)>>4;
			shift++;
		}
		*(data+length-2) = (*(data+length-1)&0x0F)<<4;
		*(data+length-1) = 0x00;
	}

	//#################### PUBLIC METHODS ####################
public:
	int find_set(int node_id) const {
		// right shift
		// & with
		// Convert to int
		num_bits_to_rightshift = (num_nodes - 1) * bits_per_node;
		return partition_vector[node_id];
	}

	void unassign_node(int node_id) {
		std::vector<int> set_seen_frequency (num_nodes, 0);

		for (int i=0; i<num_nodes; ++i) {
			set_seen_frequency[i]++;
		}
		std::vector<int>::iterator set_itr = set_seen_frequency.find(0);

		int available_set = set_itr - set_seen_frequency.begin();

		add_node_to_set(node_id, available_set);
	}

	void add_node_to_set(int node_id, int set_id) {
		// Create mask with all 1s except for 0s where node I want to change is
		// Create mask w
		for (int i=0<)
			partition_vector[node_id] = set_id;
	}

	void normalise_ids() {
		int start_num = 0;
		std::map<int, int> set_old_to_new;
		for (std::vector<int>::iterator itr = partition_vector.begin(); itr
		!= partition_vector.end(); ++itr) {

			std::map<int, int>::iterator old_set = set_old_to_new.find(*itr);
			if (old_set == set_old_to_new.end()) {
				set_old_to_new[*itr] = start_num;
				*itr = start_num;
				start_num++;
			} else {
				*itr = old_set->second;
			}
		}
	}

	int element_count() {
		return num_nodes;
	}

	int set_count() {
		std::set<int> seen_sets;

		for (int i=0; i<num_nodes; ++i) {
			int seen_set = find_set(i);
			seen_sets.insert(seen_set);
		}
		return seen_sets.size();
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
