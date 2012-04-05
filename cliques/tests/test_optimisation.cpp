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
#include <cliques/algorithms/stability.h>
#include <cliques/algorithms/space.h>
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
    lemon::SmartGraph orange_graph;
    lemon::SmartGraph::EdgeMap<double> weights(orange_graph);
    int num_samples;
    parse_arguments(ac, av, orange_graph, weights, num_samples);

    double markov_time = 1.0;
    double precision = 1e-9;
    int num_steps_per_sample = 10;

    cliques::find_full_normalised_stability func(orange_graph, weights,
            precision);
    cliques::VectorPartition initial_partition(lemon::countNodes(orange_graph));
    initial_partition.initialise_as_singletons();

    // Use a lambda to close over orange_graph so that the function passed to stochastic_climb only takes one param
    // - the partition, but has necessary access to orange_graph
    cliques::VectorPartition optimal_partition = cliques::stochastic_monotonic_climb<cliques::VectorPartition>(
            initial_partition,
            [&orange_graph] (cliques::VectorPartition partition) -> std::vector<cliques::VectorPartition> {
                std::vector<cliques::VectorPartition> alpha;
                return alpha;
//                return cliques::find_neighbours(orange_graph, partition);
            },
            cliques::Direction::ASCENT,
            [] (cliques::VectorPartition partition) -> double {
                return 1.0;
            });

    return 0;
}
