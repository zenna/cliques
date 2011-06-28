$(document).ready( function() {
	$.getJSON('js/data/ring_comm2_n6.json', function(data) {
		landscape_view = new App();
		landscape_view.setup();
		landscape = new Graph(data, landscape_view.scene, 'landscape');
		landscape.create_random_offsets(0.02);
		landscape.place_nodes();
		landscape_view.camera.target.position.x = landscape.mean_position.x;
		landscape_view.camera.target.position.y = landscape.mean_position.y;
		landscape_view.camera.target.position.z = landscape.mean_position.z;
		landscape_view.camera.position.z = 4000;

		landscape.add_edges({
			opacity:0.6
		});
		//landscape.update_edge_colours();

		graph_view = new App(300,300, 'graph_view');
		graph_view.setup();
		orig_graph = new Graph(data.graph, graph_view.scene);
		orig_graph.create_random_offsets(0);
		orig_graph.place_nodes();
		graph_view.camera.target.position.x = orig_graph.mean_position.x;
		graph_view.camera.target.position.y = orig_graph.mean_position.y;
		graph_view.camera.target.position.z = orig_graph.mean_position.z;
		graph_view.camera.position.z = 1300;
		orig_graph.add_edges({
			opacity:1.0
		});

		box = new cliques.SelectionBox(graph_view.camera.domElement);

		toolbox = new ProcessToolbox(landscape_view);
		stabilities = data.processes[0];
		basins = data.processes[1];
		toolbox.addProcess(new NodeProcess(stabilities, landscape));
		// toolbox.addProcess(new BasinProcess(basins, landscape));
		// toolbox.addProcess(new NodeProcess(louvain, landscape));
		// toolbox.addProcess(new NodeProcess(louvain2, landscape));

	}, function() {
		alert("could not load data");
	})
})