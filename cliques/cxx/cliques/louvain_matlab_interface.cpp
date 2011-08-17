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

// some redefinitions needed here as MATLAB tries to define its own symbols...
//namespace matlab {
#define CHAR16_T UINT16_T // VALID?? be careful!!
#include "mex.h"
#include "matrix.h"
//}


// GLOBAL DATA
double *data = NULL;
double precision = 1e-9; // default value here

std::vector<double> m_times;
int num_iterations = 0;

bool hierarchy = false;
int num_largest_dim = -1;

// parse function definition, needed by constructor
bool parse_arg(int nrhs, const mxArray *prhs[]) {

	//FIRST ARGUMENT: Graph
	if (nrhs > 0) {
		// number of columns should be 3 (n1,n2, weight)
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

		// get pointer to input data
		double* input_times = ((double *) mxGetPr(prhs[1]));
		// get number of columns of time vector, i.e. number of time points
		unsigned int num_times(mxGetM(prhs[1]));
		//std::cout << "number of times" << num_times << std::endl;
		//mexPrintf("Number of times %d", num_times);

		// create time vector for internal use..
		m_times.clear();
		for (int i = 0; i < num_times; ++i) {
			m_times.push_back(input_times[i]);
		}
		//std::cout << "number of times measured" << m_times.size() << std::endl;
	}

	//THIRD ARGUMENT: iterations
	if (nrhs > 2) {
		num_iterations = int((double) mxGetScalar(prhs[2]));
		std::cout << "number of iterations" << num_iterations << std::endl;
		if (num_iterations < 1) {
			return false;
		}
	}

	//FOURTH ARGUMENT: precision
	if (nrhs > 3) {
		precision = ((double) mxGetScalar(prhs[3]));
		if (precision > 1) {
			return false;
		}
	}
	// TODO: this is not implemented fully/has no effect so far.
	//FIFTH ARGUMENT: hierarchical output
	if (nrhs > 4) {
		// get buffer length and allocate buffer
		char *buf;
		mwSize buflen = mxGetN(prhs[4]) * sizeof(mxChar) + 1;
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
	if (nrhs > 5 || nrhs < 1) {
		return false;
	}
	return true;
}

// Template for reading in graph from weighted edgelist data as coming from Matlab
template<typename G, typename E>
bool read_edgelist_weighted_from_data(int num_l_dim, G &graph, E &weights) {

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
		int node1_id = data[i];
		int node2_id = data[num_l_dim + i];

		// TODO maybe there is a neater solution here
		// read in list is two-way yet undirected, but edges should only be created once
		if (node1_id > node2_id) {
			continue;
		}

		//TODO adapt for the case where unweighted graph is passed
		double weight = data[2 * num_l_dim + i];

		typename G::Edge edge = graph.addEdge(graph.nodeFromId(node1_id),
				graph.nodeFromId(node2_id));
		weights.set(edge, weight);
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
	lemon::SmartGraph::EdgeMap<double> myweights(mygraph);

	if (!read_edgelist_weighted_from_data(num_largest_dim, mygraph, myweights)) {
		mexErrMsgTxt("Error creating graph from data");
	}

	// typedef for convenience
	typedef cliques::VectorPartition partition;

	// get number of nodes and time points;
	unsigned int num_nodes = lemon::countNodes(mygraph);
	unsigned int num_times = m_times.size();

	// vector of initial partition
	partition singletons(num_nodes);
	singletons.initialise_as_singletons();

	// vector of vectors for output partitions
	std::vector<partition> optimal_partitions;

	//initialise stabilities
	std::vector<double> stability(num_iterations, 0);

	// TODO: IMPORTANT implement iterations over time correctly
	//atm this relies on the MATLAB input to be a single number
	for (int i = 0; i < num_iterations; ++i) {
		for (unsigned int j = 0; j < num_times; ++j) {

			// Logging feature
			cliques::Logging<partition> log_louvain;

			// create empty vector of partitions
			std::vector<partition> hierarchical_louvain_partitions;

			// initialise quality function
			cliques::find_linearised_normalised_stability quality(m_times[j]);

			// now run Louvain method
			stability[i] = cliques::find_optimal_partition_louvain<partition>(
					mygraph, myweights, quality,
					cliques::linearised_normalised_stability_gain(m_times[j]),
					precision, singletons, hierarchical_louvain_partitions,
					log_louvain);

			// last partition in vector == best partition
			partition best_partition = hierarchical_louvain_partitions.back();

			// store best partition
			optimal_partitions.push_back(best_partition);
		}
	}

	//****************************************************
	//----------------------------------------------------
	// Now write data back to Matlab
	//----------------------------------------------------

	// column and row dimensions
	mwSize columns = num_iterations;
	mwSize rows = num_nodes;

	/////////////////////////////////////////
	// FIRST output: stability

	// mxReal is our data-type
	plhs[0] = mxCreateDoubleMatrix(1, columns, mxREAL);
	//Get a pointer to the data space in our newly allocated memory
	double * out1 = (double*) mxGetPr(plhs[0]);

	// write out stabilities
	for (int k = 0; k < columns; ++k) {
		out1[k] = double(stability[k]);
	}

	////////////////////////////////////////
	// SECOND output: number of communities
	if (nlhs > 1) {

		plhs[1] = mxCreateDoubleMatrix(1, columns, mxREAL); //mxReal is our data-type

		//Get a pointer to the data space in our newly allocated memory
		double * out2 = (double*) mxGetPr(plhs[1]);

		// write out number of communities
		for (int k = 0; k < columns; ++k) {
			out2[k] = double(optimal_partitions[k].set_count());
		}

	}
	////////////////////////////////////////
	// THIRD output: community assignments
	if (nlhs > 2) {

		// allocate storage
		plhs[2] = mxCreateDoubleMatrix(rows, columns, mxREAL);
		double * output_tab = (double*) mxGetPr(plhs[2]);

		// write out results
		for (unsigned int iteration = 0; iteration < columns; ++iteration) {
			for (unsigned int node = 0; node < rows; ++node) {
				output_tab[iteration * num_nodes + node]
						= double(optimal_partitions[iteration].find_set(node));
			}
		}
	}

}// end namespace matlab_interface
