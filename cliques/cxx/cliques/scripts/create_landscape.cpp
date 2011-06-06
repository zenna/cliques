#include <vector>
#include <iostream>
#include <algorithm>

#include <lemon/smart_graph.h>
#include <boost/program_options.hpp>


#include <cliques/helpers.h>
#include <cliques/algorithms/all_partitions.h>
#include <cliques/algorithms/stability.h>
#include <cliques/nldr/nldr.h>
#include <cliques/structures/make_graphs.h>

#include <cliques/algorithms/space.h>
#include <cliques/algorithms/internals/internals.h>
#include <cliques/structures/vector_partition.h>

namespace po = boost::program_options;

template <typename G, typename M>
void parse_arguments(int ac, char *av[], G &graph, M &weights, int &num_samples) {
    // Declare the supported options.
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("graph,G", po::value<std::string>(),"input graph")
        ("num-samples,S", po::value<int>(),"number of samples")
    ;

    po::variables_map vm;
    po::store(po::parse_command_line(ac, av, desc), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cout << desc << "\n";
        //return 1; - exit
    }

    if (vm.count("num-samples")) {
        num_samples = vm["num-samples"].as<int>();
    }

    cliques::output("making graph");
    if (vm.count("graph")) {
        std::string filename = vm["graph"].as<std::string>();
        cliques::read_edgelist_weighted(filename, graph, weights);
    } else {
//        cliques::make_path_graph(graph, 8, weights);
//      cliques::make_ring_graph(graph, 12, weights);
      cliques::make_complete_graph(graph, 8, weights);
    }
}

int main(int ac, char* av[]) {
    typedef cliques::VectorPartition VecPartition;
    typedef boost::unordered_set<VecPartition, cliques::partition_hash,
                cliques::partition_equal> VecPartitionSet;

    lemon::SmartGraph orange_graph;
    lemon::SmartGraph::EdgeMap<double> weights(orange_graph);
    int num_samples = 100000;
    parse_arguments(ac, av, orange_graph, weights, num_samples);

    cliques::output("Finding Connected Partitions");
    cliques::NoLogging no_logging;
    cliques::Logging<VecPartition> log_all;
    VecPartitionSet all_partitions;
    cliques::find_connected_partitions(orange_graph, all_partitions, no_logging);
    cliques::output("complete size:",all_partitions.size());

    cliques::output("Creating space graph");
    lemon::SmartGraph space;
    auto map = cliques::create_space(orange_graph, all_partitions, space);
    lemon::SmartGraph::EdgeMap<float> space_weights(space);
    cliques::make_weights_from_edges(space,space_weights);

    cliques::output("Finding stabilities");
    std::vector<double> markov_times = { 1.0 };
    cliques::find_weighted_linearised_stability func(markov_times);
    std::map<int, double> stabilities;
    for (lemon::SmartGraph::NodeIt itr(space); itr != lemon::INVALID; ++itr) {
        std::vector<double> stabs;
        VecPartition p = map.right.at(itr);
        cliques::LinearisedInternals internals(orange_graph, weights, p);
        double stability = func(internals);
        stabilities[orange_graph.id(itr)] = stability;
    }
    arma::rowvec stabs_mat(stabilities.size());
    for (auto itr = stabilities.begin(); itr != stabilities.end(); ++itr) {
        stabs_mat(itr->first) = itr->second;
    }
    stabs_mat.save("stabs.mat", arma::raw_ascii);

    cliques::output("Finding distances");
    //auto X = cliques::find_geodesic_dists(space, landmark_nodes, space_weights);
    auto X = cliques::find_edit_dists(all_partitions);

    cliques::output("finding embedding");
    auto L = cliques::embed_mds(X, 3);
    arma::mat L_t = arma::trans(L);
    L_t.save("coords-full.mat", arma::raw_ascii);

    auto D_y = cliques::euclid_pairwise_dists(L_t);
    cliques::output("residual variance", cliques::residual_variance(X, D_y));

    cliques::output("number of nodes", lemon::countNodes(space));
    cliques::output("number of edges", lemon::countEdges(space));

    return 0;
}
