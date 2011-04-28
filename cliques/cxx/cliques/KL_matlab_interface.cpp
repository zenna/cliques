/*
 * matlab_interface.cpp
 *
 * compile from MATLAB in directory cxx:
 * TODO add extra options for simpler name and ensure code optimization
 *   mex  -DUSE_BOOST -I./ ./cliques/KL_matlab_interface.cpp
 *
 *  Created on: 11 Apr 2011
 *      Author: mts09
 */

#include "mex.h"
#include <lemon/smart_graph.h>
#include <map>
#include "cliques/algorithms/kernighan_lin.h"
#include "cliques/algorithms/stability.h"
#include "cliques/structures/vector_partition.h"
#include <vector>
#include <string>

double *data = NULL;
double *partition_data = NULL;
double precision = 0.000001;

// TODO enable passing time vectors..
double m_time = 1;

bool hierarchy = false;
int num_largest_dim = -1;

typedef cliques::VectorPartition partition;
partition refined_partition(1);
partition input_partition(1);

// parse function definition, needed by constructor
bool parse_arg(int nrhs, const mxArray *prhs[]) {

	//FIRST ARGUMENT: Graph
	if (nrhs > 0) {
		// number of colums should be 3 (n1,n2, weight)
		if (mxGetN(prhs[0]) != 3) {
			mexPrintf("Graph number of columns %d \n", mxGetN(prhs[0]));
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

	//THIRD ARGUMENT: partition
	if (nrhs > 2) {
		// number of rows should be 3 (n1,n2, weight)
		if (mxGetN(prhs[2]) != 1) {
			mexPrintf("Partition Number of columns %d \n", mxGetN(prhs[2]));
			return false;
		}
		// data is stored in column major ordering ordering
		partition_data = (double *) mxGetPr(prhs[2]);
	}
	//FOURTH ARGUMENT: hierarchical output
	if (nrhs > 3) {
		// get buffer length and allocate buffer
		char *buf;
		mwSize buflen = mxGetN(prhs[3]) * sizeof(mxChar) + 1;
		buf = (char*) mxMalloc(buflen);
		if (!mxGetString(prhs[3], buf, buflen)) {
			//if reading successful string is created from matlab input
			const std::string input(buf);
			const std::string comparison("h");

			// in case hierarchical output is activated print message
			if (!comparison.compare(input)) {
				hierarchy = true;
				mexPrintf("Hierarchical output from Louvain activated \n");
			}
		}
		// free read buffer
		mxFree(buf);
	}

	// SANITY CHECK for number of arguments
	if (nrhs > 4 || nrhs < 1) {
		return false;
	}
	return true;
}

// Template for reading in graph from weighted edgelist data as coming from Matlab
// TODO: adapt this to make it read "two way" files as normally used by stability code without creating double edges
// TODO: passing num_l_dim is not necessary really
template<typename G, typename E>
bool read_edgelist_weighted_from_data(double* graph_data, int num_l_dim,
		G &graph, E &weights) {

	// Find number of nodes
	// TODO: relies on correct ordering of nodes,
	// column first ordering MATLAB & node numbering starting from 0
	// but should be called from matlab anyway so maybe not important?!
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
		//TODO adapt for the case where unweighted graph is passed
		double weight = graph_data[2 * num_l_dim + i];

		// TODO maybe there is a neater solution here
		// read in list is two-way yet undirected, but edges should only be created once
		if (node1_id > node2_id) {
			continue;
		}

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

	return true;
}

template<typename G>
bool initialise_partitions(G &graph, double* part_data) {
	int number_of_nodes = lemon::countNodes(graph);
	// intialise refined partition
	refined_partition = partition(number_of_nodes);
	input_partition = partition(number_of_nodes);

	for (int i = 0; i < number_of_nodes; ++i) {
		input_partition.add_node_to_set(i, part_data[i]);
	}
	return true;
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

	if (!initialise_partitions(mygraph, partition_data)) {
		mexErrMsgTxt("Error creating partition from data");
	}

	// create time vector
	std::vector<double> markov_times;
	markov_times.push_back(m_time);

	//initialise stability
	double improvement = 0;

	// now run Louvain method
	improvement = cliques::refine_partition_kernighan_lin(mygraph, myweights,
			cliques::find_weighted_linearised_stability(markov_times),
			cliques::linearised_stability_gain_louvain(m_time),
			input_partition, refined_partition);

	//****************************************************
	//----------------------------------------------------
	// Now write data back to Matlab
	//----------------------------------------------------

	/////////////////////////////////////////
	// FIRST output: stability improvement

	// mxReal is our data-type
	plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
	//Get a pointer to the data space in our newly allocated memory
	double * out1 = (double*) mxGetPr(plhs[0]);
	out1[0] = improvement;

	////////////////////////////////////////
	// SECOND output: number of communities
	if (nlhs > 1) {

		plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL); //mxReal is our data-type

		//Get a pointer to the data space in our newly allocated memory
		double * out2 = (double*) mxGetPr(plhs[1]);

		out2[0] = double(refined_partition.set_count());

	}
	////////////////////////////////////////
	// THIRD output: community assignments
	if (nlhs > 2) {
		// get number of nodes
		int num_nodes = refined_partition.element_count();

		// allocate storage
		plhs[2] = mxCreateDoubleMatrix(num_nodes, 1, mxREAL);
		double * output_tab = (double*) mxGetPr(plhs[2]);

		// write out results
		for (unsigned int node = 0; node < num_nodes; ++node) {
			output_tab[node] = double(refined_partition.find_set(node));

		}
	}

}// end namespace matlab_interface
