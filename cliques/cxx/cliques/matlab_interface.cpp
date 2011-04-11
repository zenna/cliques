/*
 * matlab_interface.cpp
 *
 *  Created on: 11 Apr 2011
 *      Author: mts09
 */

#include "mex.h"

namespace matlab_interface {

double *data = NULL;
double precision = 0.000001;

// TODO enable passing time vectors..
double time = 1;

int display_level = -1;
int num_largest_dim = -1;

// parse function definition, needed by constructor
bool parse_arg(int nrhs, const mxArray *prhs[]) {

	//FIRST ARGUMENT: Graph
	if (nrhs > 0) {
		// number of rows should be 3 (n1,n2, weight)
		if (mxGetN(prhs[0]) != 3) {
			printf("N=%d", mxGetN(prhs[0]));
			return false;
		}
		// data is stored in column major ordering ordering
		data = (double *) mxGetPr(prhs[0]);
		// number of rows is important for indexing
		num_largest_dim = mxGetM(prhs[0]);
	}

	//SECOND ARGUMENT: time
	if (nrhs > 1) {
		time = ((double) mxGetScalar(prhs[1]));
	}

	//THIRD ARGUMENT: precision
	if (nrhs > 2) {
		if (precision > 1)
			return false;
		precision = ((double) mxGetScalar(prhs[2]));
	}
	//FOURTH ARGUMENT: hierarchical output
	if (nrhs > 3) { //TODO adapt this
		double p = (double) mxGetScalar(prhs[4]);
		if (p == 104) {
			hierarchy = true;
		} else if (p == 110) {
			hierarchy = false;
		} else {
			return false;
		}
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
		mexErrMsgTxt('Error parsing arguments');
	}
	//create new graph and weight map
	if (!read_edgelist_weighted_from_data(data, num_largest_dim, mygraph,
			myweights)) {
		mexErrMsgTxt('Error creating graph from data');
	}

	// now run louvain method
	// write data back to matlab

}

}// end namespace matlab_interface
