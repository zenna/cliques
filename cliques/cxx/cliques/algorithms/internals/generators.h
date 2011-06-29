#include <cliques/algorithms/internals/linearised_internals.h>
#include <cliques/algorithms/stability.h>

namespace cliques {

    template <typename G, typename M>
    cliques::LinearisedInternals gen(find_weighted_linearised_stability, G &graph, M &weights) {
        LinearisedInternals internals(graph, weights);
        return internals;
    }

    template <typename G, typename M, typename P>
    cliques::LinearisedInternals gen(find_weighted_linearised_stability, G &graph, M &weights, P &partition) {
        LinearisedInternals internals(graph, weights, partition);
        return internals;
    }

}
