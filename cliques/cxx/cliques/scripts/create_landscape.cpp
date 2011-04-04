#include <iostream>

#include <lemon/smart_graph.h>

#include <vector>

#include <cliques/algorithms/hsg.h>
#include <cliques/graphhelpers.h>
#include <cliques/drawing/draw.h>
#include <cliques/drawing/colour_maps.h>
#include <cliques/algorithms/partitions.h>
#include <cliques/structures/common.h>
#include <cliques/algorithms/space.h>


int main() {
    lemon::SmartGraph orange_graph;
    lemon::SmartGraph::EdgeMap<int> weights(orange_graph);

    typedef cliques::DisjointSetForest<int> DjPartition;
    cliques::read_edgelist_weighted("/home/zenna/repos/uniproj/data/triangle.edj", orange_graph, weights);

    boost::unordered_set<DjPartition,cliques::partition_hash, cliques::partition_equal> all_partitions;
    cliques::find_connected_partitions(orange_graph, all_partitions);

    lemon::SmartGraph space;
    cliques::create_space(all_partitions,space);

};
