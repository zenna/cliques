#include <iostream>

#include <lemon/smart_graph.h>

#include <vector>
#include <map>

#include <cliques/algorithms/hsg.h>
#include <cliques/graphhelpers.h>
#include <cliques/drawing/draw.h>
#include <cliques/drawing/colour_maps.h>

int main() {
    lemon::SmartGraph orange_graph;
    lemon::SmartGraph::EdgeMap<int> weights(orange_graph);

    typedef cliques::DisjointSetForest<int> DjPartition;

    cliques::read_edgelist_weighted("/home/zenna/repos/uniproj/data/celegansweighted.edj", orange_graph, weights);

    std::vector<float> markov_times;
    markov_times.push_back(0.05);
    markov_times.push_back(1.0);
    markov_times.push_back(100.0);

    lemon::SmartGraph hsg;

    std::vector<float> current_markov_time;
    current_markov_time.push_back(1.0);
    /*DjPartition best_partition = cliques::find_optimal_partition_louvain<DjPartition>
    				(orange_graph, cliques::find_linearised_stability(current_markov_time));*/

    /*cliques::find_weighted_linearised_stability<lemon::SmartGraph::EdgeMap<int> > functor(current_markov_time, weights);

    DjPartition best_partition = cliques::find_optimal_partition_louvain<DjPartition>
        				(orange_graph, functor);

    cliques::draw_graph canvas2(orange_graph);
    canvas2.add_weights(orange_graph, weights);
    canvas2.add_node_map(cliques::make_partition_colour_map<cliques::DisjointSetForest<int> >(best_partition));
    canvas2.draw("weighted.png");*/

    // Hierarchical Drawing
    cliques::draw_graph canvas(orange_graph);
    typedef cliques::DisjointSetForest<int> DjForest;
    std::vector<float> layer_stabilities;
	std::vector<DjForest> layers = cliques::create_hsg_layers<DjForest>(orange_graph, markov_times, layer_stabilities);

	int noi = 8;
	int layer_index = 0;
	for (std::vector<DjForest>::iterator itr = layers.begin(); itr != layers.end(); ++itr) {
		canvas.add_node_map(cliques::highlight_group_node_map<DjForest>(*itr,noi));
		std::ostringstream filename;
		filename << "layers" << layer_index << ".png";
		canvas.draw(filename.str());
		++layer_index;
	}

	// Hierarchical Stability Graph
	/*std::map<int,std::vector<float> > positions;
    cliques::create_hsg(orange_graph, markov_times, hsg, positions);

    cliques::draw_graph canvas(hsg);
    canvas.draw("hsg.png", "dot");*/
};
