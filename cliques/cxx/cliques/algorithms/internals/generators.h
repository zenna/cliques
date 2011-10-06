#include <cliques/algorithms/internals/linearised_internals.h>
#include <cliques/algorithms/stability.h>

namespace cliques {

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

}
