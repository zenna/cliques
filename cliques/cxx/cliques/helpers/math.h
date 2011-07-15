#pragma once

#include <iostream>
#include <vector>

extern "C" void dgpadm_(int* ideg, int* m, double* t, double* A, int* ldh,
		double* wsp, int* lwsp, int* iwsp, int* iexp, int* ns, int* iflag);

namespace cliques {

std::vector<double> exp(std::vector<double> matrix, double t, int order) {
	int ideg = 6;
	int m = order;

	int lda = order;
	int ldh = lda;
	int lwsp = 4 * ldh * ldh + ideg + 1;
	double *wsp = new double[lwsp];
	double *A = &matrix.front();
	int *iwsp = new int[ldh];

	// output
	int iexp, ns, iflag;
	dgpadm_(&ideg, &m, &t, A, &lda, wsp, &lwsp, iwsp, &iexp, &ns, &iflag);

	double *start = wsp + iexp - 1;

	std::vector<double> output;
	for (unsigned int i = 0; i < matrix.size(); ++i) {
		double a = start[i];
		output.push_back(a);
	}

	delete[] iwsp;
	delete[] wsp;

	return output;
}

std::vector<double> create_exponential_markov_times(double start_time,
		double end_time, int num_steps) {
	std::vector<double> markov_times;
	double start_log = std::log(start_time);
	double end_log = std::log(end_time);

	double increment = (end_log - start_log) / float(num_steps - 1);
	for (int i = 0; i < num_steps; ++i) {
		double current_log = start_log + i * increment;
		markov_times.push_back(std::exp(current_log));
	}
	return markov_times;
}

double discrete_gauss_kernel(int N, double T) {
	return 1.0/(N+1);
}

}
