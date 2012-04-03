#pragma once

#include <cliques/algorithms/internals/internals.h>
#include <cliques/algorithms/stability.h>
#include <cliques/algorithms/stability_corr.h>
#include <cliques/algorithms/stability_info.h>

namespace cliques {
/////////////////////////////////////////////
// NORMALISED STABILIY
/////////////////////////////////////////////
// Linearised normalised stability w/o partition
template<typename G, typename M>
cliques::LinearisedInternals gen(find_linearised_normalised_stability,
		G &graph, M &weights) {
	LinearisedInternals internals(graph, weights);
	return internals;
}

// Linearised normalised stability with partition given
template<typename G, typename M, typename P>
cliques::LinearisedInternals gen(find_linearised_normalised_stability,
		G &graph, M &weights, P &partition, P partition_unused = P(0),
		std::vector<double> null_model_vec = std::vector<double>()) {
	LinearisedInternals internals(graph, weights, partition);
	return internals;
}
//////////////////////////////////////////////
// COMBINATORIAL STABILITY
//////////////////////////////////////////////
// Linearised combinatorial stability w/o partition
template<typename G, typename M>
cliques::LinearisedInternalsComb gen(find_linearised_combinatorial_stability,
		G &graph, M &weights) {
	LinearisedInternalsComb internals(graph, weights);
	return internals;
}

// Linearised combinatorial stability with partition given
template<typename G, typename M, typename P>
cliques::LinearisedInternalsComb gen(find_linearised_combinatorial_stability,
		G &graph, M &weights, P &partition, P &partition_init, std::vector<
				double> null_model_vec= std::vector<double>()) {
	LinearisedInternalsComb
			internals(graph, weights, partition, partition_init);
	return internals;
}
/////////////////////////////////////////////
// NORMALISED CORRELATION STABILIY
/////////////////////////////////////////////
// Linearised normalised correlation stability w/o partition
template<typename G, typename M>
cliques::LinearisedInternalsCorr gen(find_linearised_normalised_corr_stability,
		G &graph, M &weights, std::vector<double> null_model_vec) {
	LinearisedInternalsCorr internals(graph, weights, null_model_vec);
	return internals;
}

// Linearised normalised correlation stability with partition given
template<typename G, typename M, typename P>
cliques::LinearisedInternalsCorr gen(find_linearised_normalised_corr_stability,
		G &graph, M &weights, P &partition, P & partition_init, std::vector<
				double> null_model_vec) {
	LinearisedInternalsCorr internals(graph, weights, partition,partition_init,null_model_vec);
	return internals;
}

/////////////////////////////////////////////
// MUTUAL INFORMATION STABILIY
/////////////////////////////////////////////
// Linearised normalised stability w/o partition
template<typename G, typename M>
cliques::LinearisedInternalsInfo gen(find_mutual_information_stability,
		G &graph, M &weights) {
	LinearisedInternalsInfo internals(graph, weights);
	return internals;
}

// Linearised normalised stability with partition given
template<typename G, typename M, typename P>
cliques::LinearisedInternalsInfo gen(find_mutual_information_stability,
		G &graph, M &weights, P &partition, P partition_unused = P(0),
		std::vector<double> null_model_vec = std::vector<double>()) {
	LinearisedInternalsInfo internals(graph, weights, partition);
	return internals;
}

}
