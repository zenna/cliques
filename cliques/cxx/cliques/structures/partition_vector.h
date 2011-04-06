/*
 * partition_vector.h
 *
 *  Created on: 6 Apr 2011
 *      Author: mts09
 */

#ifndef PARTITION_VECTOR_H_
#define PARTITION_VECTOR_H_

#include <vector>
#include <iostream>

namespace cliques {

class partition_vector {
private:
	int number_of_nodes;
	std::vector<unsigned int> partition_vector;

public:
	//#################### CONSTRUCTORS ####################

	// construct empty partition
	explicit partition_vector(int num_nodes) :
		num_nodes(num_nodes), partition_vector(num_nodes, -1) {
	}

	// construct partition from vector
	explicit partition_vector(std::vector<unsigned int> partition) :
		partition_vector(partition) {
	}
	// default destructor
	~partition_vector(){}

	bool add_node_to_set(unsigned int node_id, int comm_id){
		partition_vector[node_id] = comm_id;
	}

	int community_id(unsigned int node_id){
		return partition_vector[node_id];
	}


};

}

#endif /* PARTITION_VECTOR_H_ */
