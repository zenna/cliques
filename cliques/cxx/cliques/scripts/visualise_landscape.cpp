#include <iostream>
#include <algorithms/complexity.h>
#include <algorithms/partitions.h>
#include <helpers.h>
#include <algorithms/stability.h>
#include <algorithms/module.h>
#include <structures/disjointset.h>

#include <lemon/list_graph.h>
#include <lemon/smart_graph.h>
#include <lemon/concepts/graph_components.h>
#include <lemon/concepts/graph.h>

#include <vector>

int main() {

	// load space
	// or create space
		// Load graph from edgelist -> G
		// Create_landscape(moveset)

	// load energies
	// or Create energies
		// over markov time range

	// load
	// or generate positions (with graph viz)
	positions = get_positions();

	// map energies to colour map
	cliques::draw_graph draw(graph) //graph initialised, positioned generated;
	for time in markov time
		cliques::draw(graph, colourmap(time));

};
