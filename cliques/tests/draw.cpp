#include <iostream>
#include <lemon/smart_graph.h>
#include <lemon/concepts/graph.h>

#include <cliques/drawing/draw.h>
#include <cliques/drawing/colour_maps.h>
#include <structures/disjointset.h>
#include <cliques/helpers.h>

#include <map>

int main() {
    lemon::SmartGraph orange_graph;

    for (int i=0;i<8;++i) {
        orange_graph.addNode();
    }

    orange_graph.addEdge(orange_graph.nodeFromId(0),orange_graph.nodeFromId(1));
    orange_graph.addEdge(orange_graph.nodeFromId(0),orange_graph.nodeFromId(2));
    orange_graph.addEdge(orange_graph.nodeFromId(1),orange_graph.nodeFromId(2));
    orange_graph.addEdge(orange_graph.nodeFromId(1),orange_graph.nodeFromId(3));
    orange_graph.addEdge(orange_graph.nodeFromId(3),orange_graph.nodeFromId(4));
    orange_graph.addEdge(orange_graph.nodeFromId(3),orange_graph.nodeFromId(5));
    orange_graph.addEdge(orange_graph.nodeFromId(4),orange_graph.nodeFromId(5));

    orange_graph.addEdge(orange_graph.nodeFromId(6),orange_graph.nodeFromId(0));
    orange_graph.addEdge(orange_graph.nodeFromId(6),orange_graph.nodeFromId(2));
    orange_graph.addEdge(orange_graph.nodeFromId(6),orange_graph.nodeFromId(1));

    orange_graph.addEdge(orange_graph.nodeFromId(7),orange_graph.nodeFromId(3));
    orange_graph.addEdge(orange_graph.nodeFromId(7),orange_graph.nodeFromId(5));
    orange_graph.addEdge(orange_graph.nodeFromId(7),orange_graph.nodeFromId(4));

    //Test ColourMap
	clq::DisjointSetForest<int> partition;
    int num_elements = 8;
    for (int i =0; i<num_elements; ++i) {
    	partition.add_element(i);
    }

    partition.union_sets(0,1);
    partition.union_sets(2,3);
    partition.union_sets(3,4);
    partition.union_sets(4,5);

    clq::print_partition(partition);

    /*clq::partition_to_map<clq::DisjointSetForest<int> > p(partition);
    std::map<int,std::vector<float> > colourmap = p();

    for (std::map<int,std::vector<float> >::iterator itr = colourmap.begin();
    		itr != colourmap.end(); ++itr) {
    	std::cout << "node " << itr->first << " : ";
    	clq::print_collection(itr->second);
    }*/

    //Test drawing
    clq::draw_graph canvas(orange_graph);
    //canvas.add_node_map(clq::make_partition_colour_map<clq::DisjointSetForest<int> >(partition));
    canvas.draw("test_draw.png");
    //canvas.draw();


    return 0;
};
