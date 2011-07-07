#include <vector>
#include <iostream>

#include <algorithm>
#include <vector>
#include <map>

#include <lemon/smart_graph.h>
#include <boost/program_options.hpp>

#include <cliques/helpers.h>
#include <cliques/algorithms/all_partitions.h>
#include <cliques/algorithms/all_communities.h>
#include <cliques/algorithms/stability.h>
#include <cliques/algorithms/louvain.h>
#include <cliques/algorithms/space.h>
#include <cliques/algorithms/maxima.h>
#include <cliques/algorithms/kernighan_lin.h>
#include <cliques/algorithms/internals/internals.h>
#include <cliques/nldr/nldr.h>
#include <cliques/structures/make_graphs.h>
#include <cliques/structures/vector_partition.h>
#include <cliques/algorithms/aglob.h>

namespace po = boost::program_options;

template<typename G, typename M>
void parse_arguments(int ac, char *av[], G &graph, M &weights,
        int &num_samples, int &num_dim, std::string &filename_prefix) {
    // Declare the supported options.
    po::options_description desc("Allowed options");
    desc.add_options()("help", "produce help message")("graph,G",
            po::value<std::string>(), "input graph")("num-samples,S",
            po::value<int>(), "number of samples")("dimensions,d",
            po::value<int>(), "number of dimensions")("prefix,x",
            po::value<std::string>(), "filename prefix");

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

    if (vm.count("dimensions")) {
        num_dim = vm["dimensions"].as<int> ();
    }

    if (vm.count("prefix")) {
        filename_prefix = vm["prefix"].as<std::string> ();
    } else {
        filename_prefix = "out";
    }

    cliques::output("making graph");
    if (vm.count("graph")) {
        std::string filename = vm["graph"].as<std::string> ();
        cliques::read_edgelist_weighted(filename, graph, weights);
    } else {
        //cliques::make_path_graph(graph, 7, weights);
        //      cliques::make_ring_graph(graph, 12, weights);
        cliques::make_complete_graph(graph, 7, weights);
    }
}

int main(int ac, char* av[]) {
    typedef cliques::VectorPartition VecPartition;
    typedef boost::unordered_set<VecPartition, cliques::partition_hash,
            cliques::partition_equal> VecPartitionSet;

    lemon::SmartGraph orange_graph;
    lemon::SmartGraph space;
    lemon::SmartGraph::EdgeMap<double> space_weights(space);
    std::string filename_prefix;
    lemon::SmartGraph::EdgeMap<float> weights(orange_graph);
    int num_samples = 100000;
    int num_dim = 3;
    parse_arguments(ac, av, orange_graph, weights, num_samples, num_dim,
            filename_prefix);

    cliques::output("Find Connected Communities");
    auto communities = cliques::find_connected_communities(orange_graph);
    std::vector<double> markov_times;
    for (double t = 0.00001; t < 500.0; t = t * 1.05) {
        markov_times.push_back(t);
    }
    //    auto communities = cliques::find_community_neighbours(orange_graph, a[0]);
    cliques::output(communities.size(), "communities");

    cliques::output("Finding stabilities");
    cliques::output("time_steps", markov_times.size());
    std::ofstream stabs_file;
    stabs_file.open(filename_prefix + "_energy.mat");
    cliques::find_full_normalised_stability func(orange_graph, weights);
    std::vector<std::vector<double> > all_stabilities;
    for (unsigned int i = 0; i < markov_times.size(); ++i) {
        std::vector<double> stabilities;
        stabs_file << markov_times[i] << " ";
        for (auto itr = communities.begin(); itr != communities.end(); ++itr) {
            cliques::VectorPartition p = cliques::community_to_partition(
                    orange_graph, *itr, 0);
            double stability = func(p, 1, markov_times[i]);
            stabilities.push_back(stability);
            stabs_file << stability << " ";
        }

        all_stabilities.push_back(stabilities);
        if (i + 1 != markov_times.size()) {
            stabs_file << std::endl;
        }
    }
    stabs_file.close();

    cliques::graph_to_edgelist_file(filename_prefix + "_graph_edgelist.edj", orange_graph);

    cliques::output("Finding distances");
    auto X = cliques::find_community_edit_dists(orange_graph, communities);
    cliques::convert_dists_to_graph(space,space_weights,X, 1.0);

//    cliques::output("Finding distances Hasse");
//    cliques::create_hasse_community_space(space,communities,space);
//    std::vector<lemon::SmartGraph::Node> all_nodes;
//    for (unsigned int i = 0; i < communities.size(); ++i) {
//        all_nodes.push_back(space.nodeFromId(i));
//    }
//    cliques::make_weights_from_edges(space,space_weights);
//    auto X = cliques::convert_graph_to_geodesic_dists(space,all_nodes,space_weights);

    cliques::graph_to_edgelist_file(filename_prefix + "_landscape_edgelist.edj", space);

    cliques::output("Finding embedding");
    auto L = cliques::embed_mds(X, num_dim);
    arma::mat L_t = arma::trans(L);
    L_t.save(filename_prefix + "_coords.mat", arma::raw_ascii);

    auto D_y = cliques::euclid_pairwise_dists(L_t);
    cliques::output("residual variance", cliques::residual_variance(X, D_y));

    std::ofstream vector_file;
    vector_file.open(filename_prefix + "_partitions.mat");
    for (auto itr = communities.begin(); itr != communities.end(); ++itr) {
        cliques::VectorPartition p = cliques::community_to_partition(
                orange_graph, *itr, 0);
        int length = p.element_count();
        for (int i = 0; i < length; i++) {
            vector_file << p.find_set(i) << " ";
        }
        vector_file << std::endl;
    }

    cliques::output("Finding Probabalistic Community Basins");
    std::vector<std::map<int, std::map<int, double> >> all_basins;
    int j =0;
    for (auto stabilities = all_stabilities.begin(); stabilities != all_stabilities.end(); ++ stabilities) {
        auto basins = cliques::compute_probabalistic_basins_new(space, *stabilities);
        cliques::output("time", markov_times[j], "num_basins",basins.size());
        all_basins.push_back(basins);
        ++j;
    }
    cliques::basins_to_file(filename_prefix + "_greedy_basins.mat", all_basins, markov_times);

    cliques::output("number of nodes", lemon::countNodes(space));
    cliques::output("number of edges", lemon::countEdges(space));
    return 0;
}
