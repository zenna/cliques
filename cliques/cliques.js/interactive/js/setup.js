$(document).ready( function() {
	$.getJSON('js/data/renaud_example.js', function(data) {
		landscape_view = new App();
		landscape_view.setup();
		landscape = new Graph(data, landscape_view.scene, 'landscape');
		landscape.place_nodes();
		landscape_view.camera.target.position.x = landscape.mean_position.x;
		landscape_view.camera.target.position.y = landscape.mean_position.y;
		landscape_view.camera.target.position.z = landscape.mean_position.z;
		landscape_view.camera.position.z = 4000;

		//landscape.colour_nodes();
		landscape.add_edges();
		//landscape.update_edge_colours();

		graph_view = new App(300,300, 'graph_view');
		graph_view.setup();
		orig_graph = new Graph(data.graph, graph_view.scene);
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
		toolbox.addProcess(new NodeProcess(stabilities, landscape));
		// toolbox.addProcess(new NodeProcess(louvain, landscape));
		// toolbox.addProcess(new NodeProcess(louvain2, landscape));

	})
	function int_to_rgb(integer) {
		var rgb = [];
		var n = Math.floor(integer / 256);
		var r = integer % 256;
		var n2 = Math.floor(n / 256);
		var r2 = n % 256;
		return [n2,r2,r];
	}

	function rgb_to_int(rgb) {
		return rgb[0] * 256 * 256 + rgb[1] * 256 + rgb[2];
	}

})