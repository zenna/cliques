#pragma once

#include <cliques/quality_functions/internals/internals.h>
#include <cliques/quality_functions/stability.h>
#include <cliques/quality_functions/stability_corr.h>
#include <cliques/quality_functions/stability_info.h>

namespace clq {
/////////////////////////////////////////////
// NORMALISED STABILIY
/////////////////////////////////////////////
// Linearised normalised stability w/o partition
template<typename G, typename M>
clq::LinearisedInternals gen(find_linearised_normalised_stability,
		G &graph, M &weights) {
	LinearisedInternals internals(graph, weights);
	return internals;
}

// Linearised normalised stability with partition given
template<typename G, typename M, typename P>
clq::LinearisedInternals gen(find_linearised_normalised_stability,
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
clq::LinearisedInternalsComb gen(find_linearised_combinatorial_stability,
		G &graph, M &weights) {
	LinearisedInternalsComb internals(graph, weights);
	return internals;
}

// Linearised combinatorial stability with partition given
template<typename G, typename M, typename P>
clq::LinearisedInternalsComb gen(find_linearised_combinatorial_stability,
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
clq::LinearisedInternalsCorr gen(find_linearised_normalised_corr_stability,
		G &graph, M &weights, std::vector<double> null_model_vec) {
	LinearisedInternalsCorr internals(graph, weights, null_model_vec);
	return internals;
}

// Linearised normalised correlation stability with partition given
template<typename G, typename M, typename P>
clq::LinearisedInternalsCorr gen(find_linearised_normalised_corr_stability,
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
clq::LinearisedInternalsInfo gen(find_mutual_information_stability,
		G &graph, M &weights) {
	LinearisedInternalsInfo internals(graph, weights);
	return internals;
}

// Linearised normalised stability with partition given
template<typename G, typename M, typename P>
clq::LinearisedInternalsInfo gen(find_mutual_information_stability,
		G &graph, M &weights, P &partition, P partition_unused = P(0),
		std::vector<double> null_model_vec = std::vector<double>()) {
	LinearisedInternalsInfo internals(graph, weights, partition);
	return internals;
}

}
