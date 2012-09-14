#include <vector>
#include <iostream>

#include <algorithm>
#include <vector>
#include <map>

#include <lemon/smart_graph.h>
#include <boost/program_options.hpp>
#include <unordered_set>


#include <cliques/helpers/helpers.h>
#include <cliques/algorithms/robustness.h>
#include <cliques/helpers/make_graphs.h>
#include <cliques/structures/common.h>



namespace po = boost::program_options;

template<typename G, typename M>
void parse_arguments(int ac, char *av[], G &graph, M &weights,
        int &num_samples, int &num_dim, std::string &filename_prefix) {
    // Declare the supported options.
    po::options_description desc("Allowed options");
    desc.add_options()("help", "produce help message")
    		("graph,G",po::value<std::string>(), "input graph")
    		("num-samples,S",po::value<int>(), "number of samples")
    		("dimensions,d",po::value<int>(), "number of dimensions")
            ("prefix,x",po::value<std::string>(),"filename prefix");

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
    }
    else {
    	filename_prefix = "out";
    }

    clq::output("making graph");
    if (vm.count("graph")) {
        std::string filename = vm["graph"].as<std::string> ();
        clq::read_edgelist_weighted(filename, graph, weights);
    } else {
        //clq::make_path_graph(graph, 7, weights);
        //      clq::make_ring_graph(graph, 12, weights);
        clq::make_complete_graph(graph, 7, weights);
    }
}

int main(int ac, char* av[]) {
    typedef clq::VectorPartition VecPartition;
    typedef std::unordered_set<VecPartition, clq::partition_hash,
            clq::partition_equal> VecPartitionSet;

    lemon::SmartGraph orange_graph;
    std::string filename_prefix;
    lemon::SmartGraph::EdgeMap<double> weights(orange_graph);
    int num_samples = 100000;
    int num_dim = 3;
    parse_arguments(ac, av, orange_graph, weights, num_samples, num_dim, filename_prefix);

    clq::output("starting robustness");
    for (double t = 0.01; t < 5.0; t = t * 1.005) {
        double robustness = clq::variance_of_vi(orange_graph, weights, t, 1000);
        clq::output(t, robustness);
    }
    return 0;
}
