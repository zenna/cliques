/*
 * matlab_interface.cpp
 *
 *  Created on: 11 Apr 2011
 *      Author: mts09
 */

#include "mex.h"
#include <lemon/smart_graph.h>
#include <map>
#include "cliques/algorithms/module.h"
#include "cliques/algorithms/stability.h"
#include "cliques/structures/vector_partition.h"
#include <vector>
#include <string>

double *data = NULL;
double precision = 0.000001;

// TODO enable passing time vectors..
double m_time = 1;

bool hierarchy = false;
int num_largest_dim = -1;

// parse function definition, needed by constructor
bool parse_arg(int nrhs, const mxArray *prhs[]) {

	//FIRST ARGUMENT: Graph
	if (nrhs > 0) {
		// number of rows should be 3 (n1,n2, weight)
		if (mxGetN(prhs[0]) != 3) {
			mexPrintf("Number of columns %d", mxGetN(prhs[0]));
			return false;
		}
		// data is stored in column major ordering ordering
		data = (double *) mxGetPr(prhs[0]);
		// number of rows is important for indexing
		num_largest_dim = mxGetM(prhs[0]);
	}

	//SECOND ARGUMENT: time
	if (nrhs > 1) {
		m_time = ((double) mxGetScalar(prhs[1]));
	}

	//THIRD ARGUMENT: precision
	if (nrhs > 2) {
		if (precision > 1)
			return false;
		precision = ((double) mxGetScalar(prhs[2]));
	}
	//FOURTH ARGUMENT: hierarchical output
	if (nrhs > 3) {
		// get buffer length and allocate buffer
		char *buf;
		mwSize buflen = mxGetN(prhs[3]) * sizeof(mxChar) + 1;

		buf = (char*) mxMalloc(buflen);
		if (!mxGetString(prhs[3], buf, buflen)) {
			const std::string input(buf);
			const std::string comparison("h");
			if (!comparison.compare(input)) {
				hierarchy = true;
				mexPrintf("Hierarchical output from Louvain activated \n");
			}
		}
		mxFree(buf);
	}

	// SANITY CHECK for number of arguments
	if (nrhs > 4 || nrhs < 1) {
		return false;
	}
	return true;
}

// Template for reading in graph from weighted edgelist data as coming from Matlab
template<typename G, typename E>
bool read_edgelist_weighted_from_data(double* graph_data, int num_l_dim,
		G &graph, E &weights) {

	// Find number of nodes; relies on correct ordering of nodes,
	// column first ordering MATLAB, node numbering starting from 0
	int num_nodes = int(data[num_l_dim - 1]) + 1;

	// reserve memory space for number of nodes
	graph.reserveNode(num_nodes);

	// define Node class for convenience
	typedef typename G::Node Node;
	// mapping from id to node
	std::map<int, Node> id_to_node;

	// loop over complete list
	for (int i = 0; i < num_l_dim; ++i) {

		// get nodes and weights
		// column major ordering from MATLAB
		int node1_id = graph_data[i];
		int node2_id = graph_data[num_l_dim + i];
		float weight = graph_data[2 * num_l_dim + i];

		typename std::map<int, Node>::iterator itr = id_to_node.find(node1_id);
		Node node1, node2;

		// If the node is not in the map
		// then create node and add to map
		if (itr == id_to_node.end()) {
			node1 = graph.addNode();
			id_to_node[node1_id] = node1;
		} else {
			node1 = itr->second;
		}

		// same for node 2
		itr = id_to_node.find(node2_id);
		if (itr == id_to_node.end()) {
			node2 = graph.addNode();
			id_to_node[node2_id] = node2;
		} else {
			node2 = itr->second;
		}

		typename G::Edge edge = graph.addEdge(node1, node2);
		weights.set(edge, weight);
	}
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	// Parse arguments and return 0 if there is an error
	if (!parse_arg(nrhs, prhs)) {
		// give back zero for all outputs
		plhs[0] = mxCreateDoubleMatrix(0, 0, mxREAL);
		if (nlhs > 1)
			plhs[1] = mxCreateDoubleMatrix(0, 0, mxREAL);
		if (nlhs > 2)
			plhs[2] = mxCreateDoubleMatrix(0, 0, mxREAL);
		mexErrMsgTxt("Error parsing arguments");
	}
	//create new graph and weight map

	lemon::SmartGraph mygraph;
	lemon::SmartGraph::EdgeMap<float> myweights(mygraph);

	if (!read_edgelist_weighted_from_data(data, num_largest_dim, mygraph,
			myweights)) {
		mexErrMsgTxt("Error creating graph from data");
	}

	// typedef for convenience
	typedef cliques::VectorPartition partition;

	// create empty partition vector
	std::vector<partition> optimal_partitions;

	std::vector<double> markov_times;
	markov_times.push_back(m_time);

	double stability = 0;

	// now run Louvain method
	stability = cliques::find_optimal_partition_louvain_with_gain<partition>(
			mygraph, myweights, cliques::find_weighted_linearised_stability(
					markov_times), cliques::linearised_stability_gain_louvain(
					m_time), optimal_partitions);
	partition best_partition = optimal_partitions.back();

	// Now write data back to Matlab

	/////////////////////////////////////////
	// FIRST output: stability

	// mxReal is our data-type
	plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
	//Get a pointer to the data space in our newly allocated memory
	double * out1 = (double*) mxGetPr(plhs[0]);
	out1[0] = stability;

	////////////////////////////////////////
	// SECOND output: number of communities
	if (nlhs > 1) {

		plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL); //mxReal is our data-type

		//Get a pointer to the data space in our newly allocated memory
		double * out2 = (double*) mxGetPr(plhs[1]);

		out2[0] = double(best_partition.set_count());

	}
	////////////////////////////////////////
	// THIRD output: community assignments
	if (nlhs > 2) {
		// get number of nodes
		int num_nodes = best_partition.element_count();

		// allocate storage
		plhs[2] = mxCreateDoubleMatrix(num_nodes, 1, mxREAL);
		double * output_tab = (double*) mxGetPr(plhs[2]);

		// write out results
		for (unsigned int node = 0; node < num_nodes; ++node) {
			output_tab[node] = double(best_partition.find_set(node));

		}
	}

}// end namespace matlab_interface
