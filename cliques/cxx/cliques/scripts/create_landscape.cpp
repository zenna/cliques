#include <vector>
#include <iostream>

#include <algorithm>
#include <vector>
#include <map>

#include <lemon/smart_graph.h>
#include <boost/program_options.hpp>

#include <cliques/helpers.h>
#include <cliques/helpers/math.h>
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
        int &num_samples, int &num_dim, std::string &filename_prefix,
        bool &find_partitions, bool &create_space, bool &find_stabilities,
        bool &find_distances, bool &do_embedding, bool &find_basins,
        int &num_timesteps, double &start_time, double &end_time) {
    // Declare the supported options.
    po::options_description desc("Allowed options");
    desc.add_options()("help", "produce help message")("graph,G",
            po::value<std::string>(), "input graph")("num-samples,S",
            po::value<int>(), "number of samples")("dimensions,D",
            po::value<int>(), "number of dimensions")("prefix,x",
            po::value<std::string>(), "filename prefix")("find_partitions,p",
            "Find Partitions")("create_space,s", "Create Space")(
            "find_stabilities,r", "Find Stabilities")("find_distances,d",
            "Find Distances")("do_embedding,e", "Do Embedding")(
            "find_basins,b", "Find Basins")("time_steps", po::value<int>(),
            "Number of time steps")("start_time", po::value<double>(),
            "start time")("end_time", po::value<double>(), "end time");

    po::variables_map vm;
    po::store(po::parse_command_line(ac, av, desc), vm);
    po::notify(vm);

    find_partitions = vm.count("find_partitions") ? true : false;
    find_distances = vm.count("find_distances") ? true : false;

    if (vm.count("time_steps")) {
        num_timesteps = vm["time_steps"].as<int> ();
    }

    if (vm.count("start_time")) {
        start_time = vm["start_time"].as<double> ();
    }

    if (vm.count("end_time")) {
        end_time = vm["end_time"].as<double> ();
    }

    if (vm.count("find_distances")) {
        find_distances = true;
        find_partitions = true;
    }

    if (vm.count("create_space")) {
        create_space = true;
        find_distances = true;
        find_partitions = true;
    }

    if (vm.count("find_stabilities")) {
        find_stabilities = true;
        find_partitions = true;
    }

    if (vm.count("find_basins")) {
        find_basins = true;
        create_space = true;
        find_stabilities = true;
        find_partitions = true;
        find_distances = true;
    }

    if (vm.count("do_embedding")) {
        do_embedding = true;
        find_distances = true;
        find_partitions = true;
    }

    // Do everything if nothing specified
    if (!(find_partitions || create_space || find_stabilities || find_distances
            || do_embedding || find_basins)) {
        find_partitions = create_space = find_stabilities = find_distances
                = do_embedding = find_basins = true;
    }

    if (vm.count("help")) {
        std::cout << desc << "\n";
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

    if (vm.count("graph")) {
        cliques::output("Loading Graph");
        std::string filename = vm["graph"].as<std::string> ();
        cliques::read_edgelist_weighted(filename, graph, weights);
    } else {
        cliques::output("making default graph graph");
        //cliques::make_path_graph(graph, 7, weights);
        //      cliques::make_ring_graph(graph, 12, weights);
        cliques::make_complete_graph(graph, 8, weights);
    }
}

// TODO Need to be
// Isolate needs to remove foreign color from nodes
// Need to be able See which basin is which
// Compute basin volume
// See Node statistics

// Optimisation Algorithm
// Characterise Paths into basins
// We want an optimsiation algorithm that with a high probability will follow paths to robust maxima
// Maxima in this landscape may be structurally different, and of similar or different stability values.
// Are we going to combine our finding of maxima with robustness

// Simple algorithm
// Go to a random point in the landscape
// Make a greedy walk upwards


// Observe Basins when stable and unstable
// Size of basins
// Distribution of probabilities
// Spatial Distribution of probabilities
// Landscape hide colors
// Histogram
// Volume of the basins
// Paths into basins


// Profiling
// SmartGraph::AddEdge is the bottleneck, creates a new graph for Num Partitions * Num Timesteps, 2000000 * 500, 1 billion times, 120 billion Add Edges
// In find pbasins; create a new much bigger graph just T times, 	NUM edges * 500 = 17 billion

// Generate filename based on input as well as passed arguments
// Any base logarithms

int main(int ac, char* av[]) {
    typedef cliques::VectorPartition VecPartition;
    typedef boost::unordered_set<VecPartition, cliques::partition_hash,
            cliques::partition_equal> VecPartitionSet;

    lemon::SmartGraph orange_graph;
    std::string filename_prefix;
    lemon::SmartGraph::EdgeMap<double> weights(orange_graph);
    int num_samples = 100000, num_dim = 3, num_timesteps = 500;
    double start_time = 1e-5, end_time = 500;
    bool find_partitions, create_space, find_stabs, find_basins,
            find_distances, do_embedding;
    find_partitions = create_space = find_stabs = find_basins = find_distances
            = do_embedding = false;
    parse_arguments(ac, av, orange_graph, weights, num_samples, num_dim,
            filename_prefix, find_partitions, create_space, find_stabs,
            find_distances, do_embedding, find_basins, num_timesteps,
            start_time, end_time);

    cliques::graph_to_edgelist_file(filename_prefix + "_graph_edgelist.edj",
            orange_graph);

    VecPartitionSet all_partitions;
    if (find_partitions) {
        cliques::output("Finding Connected Partitions");
        cliques::NoLogging no_logging;
        cliques::find_connected_partitions(orange_graph, all_partitions,
                no_logging);
        cliques::output("complete size:", all_partitions.size());

        std::ofstream vector_file;
        vector_file.open(filename_prefix + "_partitions.mat");
        for (auto itr = all_partitions.begin(); itr != all_partitions.end(); ++itr) {
            int length = itr->element_count();
            for (int i = 0; i < length; i++) {
                vector_file << itr->find_set(i) << " ";
            }
            vector_file << std::endl;
        }
    }

    lemon::SmartGraph space;
    lemon::SmartGraph::EdgeMap<float> space_weights(space);
    boost::bimap<boost::bimaps::unordered_set_of<cliques::VectorPartition, cliques::partition_hash,
            cliques::partition_equal>, boost::bimaps::set_of<lemon::SmartGraph::Node> > map;
    if (create_space) {
        cliques::output("Creating space graph");
        map = cliques::create_space(orange_graph, all_partitions, space);
        //		lemon::SmartGraph::EdgeMap<float> space_weights(space);
        cliques::make_weights_from_edges(space, space_weights);
        cliques::output("number of nodes", lemon::countNodes(space));
        cliques::output("number of edges", lemon::countEdges(space));
    }

    std::vector<double> markov_times;
    std::vector<std::vector<double> > all_stabilities;
    if (find_stabs) {
        cliques::output("Finding stabilities");
        std::ofstream stabs_file;
        stabs_file.open(filename_prefix + "_energy.mat");
        cliques::output(start_time, end_time, num_timesteps);
        markov_times = cliques::create_exponential_markov_times(start_time,
                end_time, num_timesteps);
        cliques::output(markov_times.size(), "time steps");
        cliques::find_full_normalised_stability func(orange_graph, weights);
        for (unsigned int i = 0; i < markov_times.size(); ++i) {
            std::vector<double> stabilities;
            stabs_file << markov_times[i] << " ";
            for (auto itr = all_partitions.begin(); itr != all_partitions.end(); ++itr) {
                double stability = func(*itr, markov_times[i]);
                stabilities.push_back(stability);
                stabs_file << stability << " ";
            }

            all_stabilities.push_back(stabilities);
            if (i + 1 != markov_times.size()) {
                stabs_file << std::endl;
            }
        }
        stabs_file.close();
        cliques::graph_to_edgelist_file(
                filename_prefix + "_landscape_edgelist.edj", space);
    }

    arma::mat X;
    if (find_distances) {
        cliques::output("Finding distances");
        //auto X = cliques::find_geodesic_dists(space, landmark_nodes, space_weights);
        //X = cliques::find_edit_dists(all_partitions);
        X = cliques::find_split_merge_dists(all_partitions);
        X.save(filename_prefix + "_dists.mat", arma::raw_ascii);
        cliques::convert_dists_to_graph(space,space_weights,X, 1.0);
    }

    if (do_embedding) {
        cliques::output("finding embedding");
        auto L = cliques::embed_mds(X, num_dim);
        arma::mat L_t = arma::trans(L);
        L_t.save(filename_prefix + "_coords.mat", arma::raw_ascii);

        auto D_y = cliques::euclid_pairwise_dists(L_t);
        cliques::output("residual variance", cliques::residual_variance(X, D_y));
    }

    if (find_basins) {
        cliques::output("Finding Probabilistic Basins");
        std::vector<std::map<int, std::map<int, double>>>all_basins;
        int j = 0;
        for (auto stabilities = all_stabilities.begin(); stabilities
                != all_stabilities.end(); ++stabilities) {
            auto basins = cliques::compute_probabalistic_basins_new(space,
                    *stabilities);
            cliques::output("time", markov_times[j], "num_basins", basins.size());
            all_basins.push_back(basins);
            ++j;
        }

        cliques::basins_to_file(filename_prefix + "_greedy_basins.mat", all_basins,
                markov_times);
    }

    bool find_klin_basins = true;
    if (find_klin_basins) {
        cliques::output("Finding Kerninghan Lin Basins");
        std::vector<std::map<int, std::map<int, double>>>all_basins;



//        for (auto map_itr = map.left.begin(); map_itr != map.left.end(); ++map_itr) {
//            cliques::print_partition_line(*map_itr);
//        }


        for (auto time = markov_times.begin(); time != markov_times.end(); ++time) {
            lemon::SmartGraph exp_graph;
            lemon::SmartGraph::EdgeMap<double> exp_graph_weights(exp_graph);
            cliques::graph_to_exponential_graph(orange_graph, weights,exp_graph, exp_graph_weights, *time);

            auto basins = cliques::compute_kernighan_lin_basins(orange_graph,weights,
                    cliques::find_linearised_normalised_stability(*time),
                    cliques::linearised_normalised_stability_gain(*time),
                    map,
                    space,
                    *time,
                    all_partitions);
            all_basins.push_back(basins);
        }
        cliques::basins_to_file(filename_prefix + "_klin_basins.mat", all_basins,
                markov_times);
    }
    return 0;
}

