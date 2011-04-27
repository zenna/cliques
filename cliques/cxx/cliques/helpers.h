#ifndef CLIQUES_HELPERS_H
#define CLIQUES_HELPERS_H

#include <iostream>
#include <cliques/structures/common.h>
#include <map>

#define N 128
#define MAX_LINE_LENGTH 1000
#define SET_MAX_SIZE 63

#define NUM_NODES 16
#define BITS_PER_NODE 4

namespace cliques {
void reorder_map(std::map<int,int> &partition_map) {
	int start_part = 0;
	std::map<int,int> temp_partition_map = partition_map;

	for (std::map<int,int>::iterator itr = temp_partition_map.begin(); itr != temp_partition_map.end(); ++itr) {
		int reference_part = itr->second;
		std::vector<int> to_remove;

		for (std::map<int,int>::iterator new_itr = temp_partition_map.begin(); new_itr != temp_partition_map.end(); ++new_itr) {
			if (new_itr->second == reference_part) {
				partition_map[new_itr->first] = start_part;

				if (itr->first != new_itr->first)
				   to_remove.push_back(new_itr->first);
			}
		}
		for (unsigned int i =0; i< to_remove.size(); ++i) {
			temp_partition_map.erase(to_remove[i]);
		}
		start_part++;
	}
}

//Sets the value of a bitset to value at offset
void set_bitset(std::bitset<N> &bit_array, int &offset, int value) {
	std::bitset<BITS_PER_NODE> temp_bit_array = value;
	for (int i =0;i<BITS_PER_NODE;i++) {
		int j = offset + i;
		bit_array[j] = temp_bit_array[i];
	}
	offset = offset + BITS_PER_NODE;
}

int get_bitset(std::bitset<N> &bit_array, int &offset) {
	std::bitset<BITS_PER_NODE> temp_bit_array;
	for (int i =0;i<BITS_PER_NODE;i++) {
		int j = offset + i;
		temp_bit_array[i] = bit_array[j];
	}
	offset = offset + BITS_PER_NODE;

	return int(temp_bit_array.to_ulong());
}

// CONVERSIONS --------------------------------------------

std::map<int,int> bitset_to_map(std::bitset<N> &bit_array) {
	int offset = 0;
	std::map<int,int> partition_map;

	for (int i = 0; i < NUM_NODES; i++) {
		int node = get_bitset(bit_array,offset);
		int part = get_bitset(bit_array,offset);
		partition_map[node] = part;
	}

	return partition_map;
};

std::bitset<N> map_to_bitset(std::map<int,int> &partition_map) {
	int offset = 0;
	std::bitset<N> bit_array;

	for (std::map<int, int>::iterator itr = partition_map.begin(); itr != partition_map.end(); ++itr) {
		set_bitset(bit_array, offset, itr->first);
		set_bitset(bit_array, offset, itr->second);
	}

	return bit_array;
}

std::vector<std::vector<int> > bitset_to_vector(std::bitset<N> &bit_array) {
	int offset = 0;
	std::vector<std::vector<int > > p;

	int first_part = 0;
	std::vector<int> temp_vec;
	for (int i = 0; i< NUM_NODES; i++) {
		int part_index = get_bitset(bit_array,offset);
		int node_index = get_bitset(bit_array,offset);
		if (part_index != first_part) {
			p.push_back(temp_vec);
			temp_vec.clear();
			first_part++;
		}
		temp_vec.push_back(node_index);
	}
	p.push_back(temp_vec);

	return p;
};

std::vector<std::vector<int> > map_to_vector(std::map<int,int> &partition_map) {
	int max_part_seen = -1;
	std::vector<std::vector<int> > temp_partition_vector;
	for (std::map<int,int>::iterator itr = partition_map.begin(); itr != partition_map.end(); ++itr) {

		int part_num_diff = itr->second - max_part_seen;
		if (part_num_diff> 0) {
			for (int i=0;i<part_num_diff;++i) {
				std::vector<int> a;
				temp_partition_vector.push_back(a);
			}
			max_part_seen = itr->second;
		}
		temp_partition_vector[itr->second].push_back(itr->first);
	}
return temp_partition_vector;
}

/*std::bitset<N> vector_to_bitset(vec_2d_ints &p) {

	std::bitset<N> temp_bitset;
	int offset = 0;
	for (unsigned int i = 0;i<p.size();++i) {
		for (unsigned int j = 0; j < p[i].size(); ++j) {
			set_bitset(temp_bitset,offset,i);
			set_bitset(temp_bitset,offset,p[i][j]);
		}
	}

	return temp_bitset;
};*/

std::set<std::set<int> > bitset_to_set(std::bitset<N> &bit_array) {
	int offset = 0;
	std::set<std::set<int > > p;

	int first_part = 0;
	std::set<int> temp_vec;
	for (int i = 0; i<16; i++) {
		int part_index = get_bitset(bit_array,offset);
		int node_index = get_bitset(bit_array,offset);
		if (part_index != first_part) {
			p.insert(temp_vec);
			temp_vec.clear();
			first_part++;
		}
		temp_vec.insert(node_index);
	}
	p.insert(temp_vec);

	return p;
}

std::vector<char> string_to_char_p(std::string str) {
	std::vector<char> writable(str.size() + 1);
	std::copy(str.begin(), str.end(), writable.begin());
	return writable;
}

// PRINTING --------------------------------------------------

void print_map(std::map<int,int> my_map, bool add_one = false) {
	for (std::map<int,int>::iterator itr = my_map.begin(); itr != my_map.end(); ++itr) {
		int node = itr->first;
		int part = itr->second;
		if (add_one) {
			node++;
			part++;
		}
		std::cout << node << ":" << part << " ";
	}
	std::cout << "\n";
}

// Prints  umap to screen
void print_umap(cliques::umap &my_umap) {
	for (cliques::umap::iterator itr = my_umap.begin(); itr != my_umap.end(); ++itr) {
		std::cout << itr->second.key << std::endl;
	}
}

void print_2d_vector (std::vector<std::vector<int> > my_vector) {
	std::vector<std::vector<int> >::iterator itr;
	for ( itr = my_vector.begin(); itr != my_vector.end(); ++itr ) {
		std::vector<int>::iterator new_itr;
		for ( new_itr = itr->begin(); new_itr != itr->end(); ++new_itr ) {
			std::cout << (*new_itr) << " ";
		}
		std::cout << "\n";
	}
}

template <class T>
void print_collection (T collection) {
	for (class T::iterator itr = collection.begin(); itr != collection.end(); ++itr) {
		std::cout << *itr << ", ";
	}
	std::cout << std::endl;
}

int find_degree(lemon::ListGraph &graph, int &node_id) {
	int count = 0;
	lemon::ListGraph::Node n = graph.nodeFromId(node_id);

	for(lemon::ListGraph::IncEdgeIt e(graph, n); e!= lemon::INVALID; ++e) {
		++count;
	}
	return count;
}

float find_weight(lemon::ListGraph &graph, int node1_id, int node2_id) {
	lemon::ListGraph::Node n1 = graph.nodeFromId(node1_id);
	lemon::ListGraph::Node n2 = graph.nodeFromId(node2_id);

	for(lemon::ListGraph::IncEdgeIt e(graph, n1); e!= lemon::INVALID; ++e) {
		std::cout << "looking at node1 " << graph.id(n1) << std::endl;
		std::cout << "looking at node2 " << graph.id(n2) << std::endl;
		std::cout << "running: " << graph.id(graph.runningNode(e)) << std::endl;
		std::cout << "base: " << graph.id(graph.baseNode(e)) << std::endl;
		if (graph.runningNode(e) == n2) {

			std::cout << "one" << std::endl;
			return 1.0;
		}
	}
	std::cout << "hero" << std::endl;
	return 0.0;

}

template <typename P>
void print_partition(P &partition) {
	for (typename P::PartIterator pitr = partition.begin(); pitr != partition.end(); ++pitr) {
		std::cout << "[";
		for (typename P::NodeIterator nitr = pitr.begin(); nitr != pitr.end(); ++nitr ) {
			std::cout << " " << *nitr << " ";
		}
		std::cout << "]" << std::endl;
	}
}

template <typename P>
void print_partition_list(P &partition) {
    int length = partition.element_count();
    for (int i = 0; i < length; i++) {
        std::cout << i << "->" << partition.find_set(i) << std::endl;
    }
}

// Type conversion
template <class T>
inline std::string to_string (const T& t)
{
	std::stringstream ss;
	ss << t;
	return ss.str();
}

template <class T>
inline std::string to_char_pointer (const T& t)
{
	std::stringstream ss;
	ss << t;
	return ss.str();
}

} //namespace cliques

#endif //CLIQUES_STRUCTURES_H
