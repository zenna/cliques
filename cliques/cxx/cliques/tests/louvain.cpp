#include <iostream>
#include <algorithms/complexity.h>
#include <algorithms/partitions.h>
#include <helpers.h>
#include <algorithms/stability.h>
#include <algorithms/modularity.h>

#include <drawing/draw.h>
#include <drawing/colour_maps.h>

#include <algorithms/module.h>
#include <structures/disjointset.h>

#include <lemon/list_graph.h>
#include <lemon/smart_graph.h>
#include <lemon/concepts/graph_components.h>
#include <lemon/concepts/graph.h>
#include <lemon/connectivity.h>


#include <vector>

#include <graphhelpers.h>


//TODO
//REMOVE NEIGHBOURS TO SELF
int main() {
    typedef cliques::DisjointSetForest<int> DjPartition;
    lemon::SmartGraph orange_graph;
    lemon::SmartGraph::EdgeMap<int> weights(orange_graph);

    cliques::read_edgelist_weighted("/home/zenna/repos/uniproj/data/celegansweighted.edj", orange_graph, weights);

    std::vector<float> current_markov_time;
    current_markov_time.push_back(1.0);

    std::vector<DjPartition> optimal_partitions;
    DjPartition best_partition = cliques::find_optimal_partition_louvain<DjPartition>
    	(orange_graph, weights, cliques::find_linearised_stability(current_markov_time), optimal_partitions);

    cliques::print_partition(optimal_partitions.back());

    //Drawing
    float start = -10;
    std::vector<float> energies;
    for (lemon::SmartGraph::EdgeIt e(orange_graph); e != lemon::INVALID; ++e) {
    	energies.push_back(start);
    	start += 12.0;
    	std::cout << orange_graph.id(e) << std::endl;
    }

    cliques::draw_graph canvas(orange_graph);
    canvas.add_node_map(cliques::make_partition_colour_map<cliques::DisjointSetForest<int> >(best_partition));
    canvas.add_edge_map(cliques::make_energy_edge_colour_map(energies));
    canvas.draw("test_louvain_out");
    return 0;
};
