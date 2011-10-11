#include <cliques/algorithms/internals/linearised_internals.h>
#include <cliques/algorithms/internals/linearised_internals_comb.h>
#include <cliques/algorithms/internals/linearised_internals_corr.h>
#include <cliques/algorithms/stability.h>

namespace cliques {
	/////////////////////////////////////////////
	// NORMALISED STABILIY
	/////////////////////////////////////////////
	// Linearised normalised stability w/o partition
    template <typename G, typename M>
    cliques::LinearisedInternals gen(find_linearised_normalised_stability, G &graph, M &weights) {
        LinearisedInternals internals(graph, weights);
        return internals;
    }

    // Linearised normalised stability with partition given
    template <typename G, typename M, typename P>
    cliques::LinearisedInternals gen(find_linearised_normalised_stability, G &graph, M &weights, P &partition, P &partition_unused) {
        LinearisedInternals internals(graph, weights, partition);
        return internals;
    }
    //////////////////////////////////////////////
    // COMBINATORIAL STABILITY
    //////////////////////////////////////////////
    // Linearised combinatorial stability w/o partition
    template <typename G, typename M>
    cliques::LinearisedInternalsComb gen(find_linearised_combinatorial_stability, G &graph, M &weights) {
        LinearisedInternalsComb internals(graph, weights);
        return internals;
    }

    // Linearised combinatorial stability with partition given
    template <typename G, typename M, typename P>
    cliques::LinearisedInternalsComb gen(find_linearised_combinatorial_stability, G &graph, M &weights, P &partition, P &partition_init) {
        LinearisedInternalsComb internals(graph, weights, partition, partition_init);
        return internals;
    }
	/////////////////////////////////////////////
	// NORMALISED CORRELATION STABILIY
	/////////////////////////////////////////////
	// Linearised normalised correlation stability w/o partition
    template <typename G, typename M>
    cliques::LinearisedInternalsCorr gen(find_linearised_normalised_corr_stability, G &graph, M &weights) {
        LinearisedInternalsCorr internals(graph, weights);
        return internals;
    }

    // Linearised normalised correlation stability with partition given
    template <typename G, typename M, typename P>
    cliques::LinearisedInternalsCorr gen(find_linearised_normalised_corr_stability, G &graph, M &weights, P &partition, P &partition_unused) {
        LinearisedInternalsCorr internals(graph, weights, partition);
        return internals;
    }


}
