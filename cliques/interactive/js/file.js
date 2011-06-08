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
	
	$.getJSON('js/data/example.js', function(data) {
		// need to be able to
		// change coordinates
		// render in just two coordinates
		// animate
		
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
		
		// 
		// camera.position.z = 4000;
		// camera.target.position.x = mean_position.x;
		// camera.target.position.y = mean_position.y;
		// camera.target.position.z = mean_position.z;
		
		// landscape.colour_nodes();
		// landscape.add_edges();
		// landscape.add_to_scene(app.scene);
// 		
  		// app.setup(data.coords);
  		// app.add_edges(data.edges);
  		// app.add_colours(data.energies);
	})
	
})