#include <vector>
#include <iostream>
#include <algorithm>
#include <unordered_map>
#include <functional>

#include <lemon/smart_graph.h>
#include <boost/program_options.hpp>

#include <cliques/helpers/helpers.h>
#include <cliques/helpers/graphhelpers.h>
#include <cliques/helpers/make_graphs.h>
#include <cliques/structures/vector_partition.h>
#include <cliques/quality_functions/stability.h>
#include <cliques/landscapes/landscape_space.h>
#include <cliques/algorithms/hill_climbing.h>

namespace po = boost::program_options;

template<typename G, typename M>
void parse_arguments(int ac, char *av[], G &graph, M &weights, int &num_samples) {
    // Declare the supported options.
    po::options_description desc("Allowed options");
    desc.add_options()("help", "produce help message")("graph,G",
            po::value<std::string>(), "input graph")("num-samples,S",
            po::value<int>(), "number of samples");

    po::variables_map vm;
    po::store(po::parse_command_line(ac, av, desc), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cout << desc << "\n";
        //return 1; - exit
    }

    if (vm.count("num-samples")) {
        num_samples = vm["num-samples"].as<int> ();
    }

    cliques::output("making graph");
    if (vm.count("graph")) {
        std::string filename = vm["graph"].as<std::string> ();
        cliques::read_edgelist_weighted(filename, graph, weights);
    } else {
        //      cliques::make_complete_graph(graph, 8, weights);
        //      cliques::make_ring_graph(graph, 12, weights);
        cliques::make_complete_graph(graph, 6, weights);
    }
}

int main(int ac, char* av[]) {
    typedef typename cliques::VectorPartition vecPart;
    typedef typename std::unordered_set<vecPart, cliques::partition_hash, cliques::partition_equal> partition_set;
    lemon::SmartGraph orange_graph;
    lemon::SmartGraph::EdgeMap<double> weights(orange_graph);
    int num_samples;
    parse_arguments(ac, av, orange_graph, weights, num_samples);

    double markov_time = 1.0;
    double precision = 1e-9;

    // std::vector<int> alpha = [0,1,2,0,0,1];

    vecPart initial_partition(lemon::countNodes(orange_graph));
    initial_partition.initialise_as_singletons();

    cliques::find_linearised_normalised_stability quality_func(markov_time);
    // cliques::find_full_normalised_stability quality_func(orange_graph, weights, precision);

    // Use a lambda to close over orange_graph so that the function passed to stochastic_climb only takes one param
    // - the partition, but has necessary access to orange_graph
    std::mt19937 prng_engine;

    cliques::output("starting to roll");
    for (int i=0;i<1000;++i) {
    vecPart optimal_partition = cliques::stochastic_monotonic_climb
            <vecPart, std::vector<vecPart>>
            (initial_partition,
            [&orange_graph] (vecPart partition) -> std::vector<vecPart> {
                return cliques::find_neighbours2(orange_graph, partition);
            },
            cliques::Direction::ASCENT,
            [&quality_func, &orange_graph, &weights] (vecPart partition) -> double {
                return quality_func(orange_graph, partition, weights);
            },
            prng_engine);

    cliques::print_partition_line(optimal_partition);

    }

    return 0;
}
