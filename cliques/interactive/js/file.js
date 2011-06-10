$(document).ready( function() {
	function handleFileSelect(evt) {
		var files = evt.target.files; // FileList object
		f = files[0];
		if (f) {
			file_to_array(f, parseFloat, app.setup);
		} else {
			alert("no file");
		}
	}

	function handle_energy_files(evt) {
		var files = evt.target.files; // FileList object
		f = files[0];
		if (f) {
			file_to_array(f, parseFloat, app.add_colours);
		} else {
			alert("no file");
		}
	}

	function handle_edges_files(evt) {
		var files = evt.target.files; // FileList object
		f = files[0];
		if (f) {
			file_to_array(f, parseInt, app.add_edges);
		} else {
			alert("no file");
		}
	}

	function file_to_array(f, parser, callback) {
		var edges_vectors = [];
		var reader = new FileReader();
		reader.onload = function(e) {
			coords = e.target.result.split("\n");
			for (var i=0;i<coords.length;++i) {
				var vector = coords[i].split(" ");
				var float_vector = [];
				for (var j=0;j<vector.length;++j) {
					if (vector[j] != " " && vector[j] != "") {
						float_vector.push(parser(vector[j]));
					}
				}
				if (float_vector.length) {
					edges_vectors.push(float_vector);
				}
			}
			callback.call(app, edges_vectors);
		}
		reader.readAsText(f);
		return edges_vectors;
	}

	document.getElementById('files').addEventListener('change',
	handleFileSelect, false);

	document.getElementById('energy_files').addEventListener('change',
	handle_energy_files, false);

	document.getElementById('edges_files').addEventListener('change',
	handle_edges_files, false);
	
	$.getJSON('js/data/complete_example.js', function(data) {
		// need to be able to
		// change coordinates
		// render in just two coordinates
		// animate
		app = new App();
		app.setup();
// 		
		landscape = new Graph(data, app.scene);
		landscape.place_nodes();
		app.camera.target.position.x = landscape.mean_position.x;
		app.camera.target.position.y = landscape.mean_position.y;
		app.camera.target.position.z = landscape.mean_position.z;
		app.camera.position.z = 4000;
		
		landscape.colour_nodes();
		landscape.add_edges();
		landscape.update_edge_colours();
	})
	
	$.getJSON('js/data/complete_example.js', function(data) {
		// need to be able to
		// change coordinates
		// render in just two coordinates
		// animate
		
		// app2 = new App();
		// app2.setup();
// // 		
		// landscape2 = new Graph(data, app2.scene);
		// landscape2.place_nodes();
		// app2.camera.target.position.x = landscape2.mean_position.x;
		// app2.camera.target.position.y = landscape2.mean_position.y;
		// app2.camera.target.position.z = landscape2.mean_position.z;
		// app2.camera.position.z = 4000;
// // 		
		// landscape2.colour_nodes();
		// landscape2.add_edges();
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