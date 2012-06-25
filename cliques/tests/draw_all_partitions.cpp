#include <iostream>

#include <lemon/smart_graph.h>

#include <vector>

#include <cliques/algorithms/hsg.h>
#include <cliques/graphhelpers.h>
#include <cliques/drawing/draw.h>
#include <cliques/drawing/colour_maps.h>
#include <cliques/algorithms/partitions.h>
#include <cliques/structures/common.h>

int main() {
    lemon::SmartGraph orange_graph;
    lemon::SmartGraph::EdgeMap<int> weights(orange_graph);

    typedef clq::DisjointSetForest<int> DjPartition;
    clq::read_edgelist_weighted("/home/zenna/repos/uniproj/data/triangle.edj", orange_graph, weights);

    std::unordered_set<DjPartition,clq::partition_hash, clq::partition_equal> all_partitions;

    clq::find_connected_partitions(orange_graph, all_partitions);
    clq::draw_graph canvas(orange_graph);

    int index = 0;
    for (std::unordered_set<DjPartition,clq::partition_hash, clq::partition_equal>::iterator itr = all_partitions.begin();
    		itr != all_partitions.end(); ++itr) {
    	clq::DisjointSetForest<int> wow = *itr;
    	canvas.add_node_map(clq::make_partition_colour_map<DjPartition>(wow));
    	std::ostringstream filename;
    	filename << "partition" << index << ".png";
    	canvas.draw(filename.str());
    	++index;
    }
};
