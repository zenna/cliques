/*
 * matlab_interface.cpp
 *
 * compile from MATLAB in directory cxx:
 * TODO add extra options for simpler name and ensure code optimization
 *   mex -DUSE_BOOST -I./ ./cliques/louvain_matlab_interface.cpp
 *
 *  Created on: 11 Apr 2011
 *      Author: mts09
 */

#include <lemon/smart_graph.h>
#include <map>
#include <cliques/algorithms/louvain.h>
#include <cliques/algorithms/stability.h>
#include <cliques/structures/vector_partition.h>
#include <vector>
#include <string>
#include <stdlib.h>

namespace matlab {
#define CHAR16_T UINT16_T // VALID?? be careful!!
#include "mex.h"
#include "matrix.h"
}






double *data = NULL;
double precision = 1e-9;

// TODO enable passing time vectors..
double m_time = 1;

bool hierarchy = false;
int num_largest_dim = -1;

// parse function definition, needed by constructor
bool parse_arg(int nrhs, const matlab::mxArray *prhs[]) {

	//FIRST ARGUMENT: Graph
	if (nrhs > 0) {
		// number of columns should be 3 (n1,n2, weight)
		if (mxGetN(prhs[0]) != 3) {
			matlab::mexPrintf("Number of columns %d", matlab::mxGetN(prhs[0]));
			return false;
		}
		// data is stored in column major ordering ordering
		data = (double *) matlab::mxGetPr(prhs[0]);
		// number of rows is important for indexing
		num_largest_dim = matlab::mxGetM(prhs[0]);
	}

	//SECOND ARGUMENT: time
	if (nrhs > 1) {
		m_time = ((double) matlab::mxGetScalar(prhs[1]));
	}

	//THIRD ARGUMENT: precision
	if (nrhs > 2) {
		if (precision > 1)
			return false;
		precision = ((double) matlab::mxGetScalar(prhs[2]));
	}
	//FOURTH ARGUMENT: hierarchical output
	if (nrhs > 3) {
		// get buffer length and allocate buffer
		char *buf;
		matlab::mwSize buflen = matlab::mxGetN(prhs[3]) * sizeof(matlab::mxChar) + 1;
		buf = (char*) matlab::mxMalloc(buflen);
		if (!matlab::mxGetString(prhs[3], buf, buflen)) {
			//if reading successful string is created from matlab input
			const std::string input(buf);
			const std::string comparison("h");

			// in case hierarchical output is activated print message
			if (!comparison.compare(input)) {
				hierarchy = true;
				matlab::mexPrintf("Hierarchical output from Louvain activated \n");
			}
		}
		// free read buffer
		matlab::mxFree(buf);
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

	// Find number of nodes
	// TODO: relies on correct ordering of nodes,
	// column first ordering MATLAB & node numbering starting from 0
	// but should be called from matlab anyway so maybe not important?!
	int num_nodes = int(data[num_l_dim - 1]) + 1;

	// reserve memory space for number of nodes
	graph.reserveNode(num_nodes);
	// now node id internal <> external
	for (int i = 0; i < num_nodes; ++i) {
		graph.addNode();
	}

	// define Node class for convenience
	typedef typename G::Node Node;

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

		typename G::Edge edge = graph.addEdge(graph.nodeFromId(node1_id),
				graph.nodeFromId(node2_id));
		weights.set(edge, weight);
	}

	return true;
}

void mexFunction(int nlhs, matlab::mxArray *plhs[], int nrhs, const matlab::mxArray *prhs[]) {

	// Parse arguments and return 0 if there is an error
	if (!parse_arg(nrhs, prhs)) {
		// give back zero for all outputs
		plhs[0] = matlab::mxCreateDoubleMatrix(0, 0, matlab::mxREAL);
		if (nlhs > 1)
			plhs[1] = matlab::mxCreateDoubleMatrix(0, 0, matlab::mxREAL);
		if (nlhs > 2)
			plhs[2] = matlab::mxCreateDoubleMatrix(0, 0, matlab::mxREAL);
		matlab::mexErrMsgTxt("Error parsing arguments");
	}

	//create new graph and weight map
	lemon::SmartGraph mygraph;
	lemon::SmartGraph::EdgeMap<float> myweights(mygraph);

	if (!read_edgelist_weighted_from_data(data, num_largest_dim, mygraph,
			myweights)) {
		matlab::mexErrMsgTxt("Error creating graph from data");
	}

	// typedef for convenience
	typedef cliques::VectorPartition partition;

    // vector of initial partition
    partition singletons(lemon::countNodes(mygraph));
    singletons.initialise_as_singletons();
    
    // Logging feature
    cliques::Logging<partition> log_louvain;
    
	// create empty vector of partitions
	std::vector<partition> optimal_partitions;

	//initialise stability
	double stability = 0;
    cliques::find_linearised_normalised_stability quality(m_time);

    // now run Louvain method
    stability = cliques::find_optimal_partition_louvain<partition>(mygraph,
            myweights, quality,
            cliques::linearised_normalised_stability_gain(m_time), precision, singletons,
            optimal_partitions, log_louvain);

	// last partition in vector == best partition
	partition best_partition = optimal_partitions.back();

	//****************************************************
	//----------------------------------------------------
	// Now write data back to Matlab
	//----------------------------------------------------

	/////////////////////////////////////////
	// FIRST output: stability

	// mxReal is our data-type
	plhs[0] = matlab::mxCreateDoubleMatrix(1, 1, matlab::mxREAL);
	//Get a pointer to the data space in our newly allocated memory
	double * out1 = (double*) matlab::mxGetPr(plhs[0]);
	out1[0] = stability;

	////////////////////////////////////////
	// SECOND output: number of communities
	if (nlhs > 1) {

		plhs[1] = matlab::mxCreateDoubleMatrix(1, 1, matlab::mxREAL); //mxReal is our data-type

		//Get a pointer to the data space in our newly allocated memory
		double * out2 = (double*) matlab::mxGetPr(plhs[1]);

		out2[0] = double(best_partition.set_count());

	}
	////////////////////////////////////////
	// THIRD output: community assignments
	if (nlhs > 2) {
		// get number of nodes
		int num_nodes = best_partition.element_count();

		// allocate storage
		plhs[2] = matlab::mxCreateDoubleMatrix(num_nodes, 1, matlab::mxREAL);
		double * output_tab = (double*) matlab::mxGetPr(plhs[2]);

		// write out results
		for (unsigned int node = 0; node < num_nodes; ++node) {
			output_tab[node] = double(best_partition.find_set(node));

		}
	}

}// end namespace matlab_interface
