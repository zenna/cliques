function Graph(data, scene) {
	this.partitions = data.partitions;
	this.coords = data.coords;
	this.energies = data.energies;
	this.edges = data.edges;
	this.mean_position = {
		x:0.0,
		y:0.0,
		z:0.0
	};
	this.scene = scene;
	this.materials = [], this.nodes = []
	this.num_dim = this.coords[0].length;

	this.place_nodes = function(dim1, dim2, dim3) {
		var materials = this.materials;
		var nodes = this.nodes;
		var coords = this.coords;
		var mean_position = this.mean_position;

		var geometry = new THREE.Geometry();
		var colours = [];

		for ( var i = 0; i < coords.length; ++i) {
			x = 2000 * Math.random() - 1000;
			y = 2000 * Math.random() - 1000;
			z = 2000 * Math.random() - 1000;
			vector = new THREE.Vector3( x, y, z );
			geometry.vertices.push( new THREE.Vertex( vector ) );
			colours[i] = new THREE.Color( 0x0084BD );
			var hue  = 0.1 * (i % 9);
			colours[ i ].setHSV(hue, 1.0, 1.0 );
		}

		this.node_id_colours = [];
		for ( var i = 0; i < coords.length; ++i) {
			var colour = this.int_to_rgb(i);
			var hex_r = colour[0].toString(16);
			var hex_g = colour[1].toString(16);
			var hex_b = colour[2].toString(16);

			var colour_hex = hex_r + hex_g + hex_b;
			this.node_id_colours[i] = new THREE.Color( parseInt(colour_hex, 16));
		}

		geometry.colors = colours;

		var material = new THREE.ParticleBasicMaterial({
			size: 40	,
			vertexColors: true
		} );
		material.color.setHSV( 0.0, 0.0, 1.0 );

		var particles = new THREE.ParticleSystem( geometry, material );
		particles.sortParticles = true;
		particles.updateMatrix();
		this.nodes = geometry;
		this.particles = particles;

		this.move_nodes();
		this.add_to_scene([particles]);
	}
	this.int_to_rgb = function(integer) {
		var n = Math.floor(integer / 256);
		var r = integer % 256;
		var n2 = Math.floor(n / 256);
		var r2 = n % 256;
		return [n2,r2,r];
	}
	this.rgb_to_int = function(rgb) {
		return rgb[0] * 256 * 256 + rgb[1] * 256 + rgb[2];
	}
	this.move_nodes = function( dim1, dim2, dim3) {

		var dim1 = typeof dim1  == "undefined" ? 0 : dim1;
		var dim2 = typeof dim2  == "undefined" ? 1 : dim2;
		var dim3;
		if (this.num_dim > 2) {
			dim3 = typeof dim3  == "undefined" ? 2 : dim3;
		} else {
			dim3 = undefined;
		}

		var materials = this.materials;
		var nodes = this.nodes;
		var coords = this.coords;
		var mean_position = this.mean_position;

		for ( var i = 0; i < nodes.vertices.length; ++i) {
			var x = coords[i][dim1] * 2000 - 1000;
			var y = coords[i][dim2] * 2000 - 1000;
			var z;
			if (typeof dim3  == "undefined") {
				z = coords[0][dim1] * 2000 - 1000;
			} else {
				z = coords[i][dim3] * 2000 - 1000;
			}
			nodes.vertices[i].position.x = x;
			nodes.vertices[i].position.y = y;
			nodes.vertices[i].position.z = z;
			nodes.__dirtyVertices = true;

			mean_position.x = mean_position.x + nodes.vertices[i].position.x
			mean_position.y = mean_position.y + nodes.vertices[i].position.y;
			mean_position.z = mean_position.z + nodes.vertices[i].position.z;
		}

		mean_position.x = mean_position.x / coords.length;
		mean_position.y = mean_position.y / coords.length;
		mean_position.z = mean_position.z / coords.length;
		if (this.lines) {
			this.move_edges();
		}
	}
	this.colour_nodes = function() {
		var energy_vectors = this.energies;
		var min_value = 1e9;
		var max_value = -1e9;
		// Find range
		for (i = 0; i < energy_vectors[0].length; i++) {
			var energy = energy_vectors[0][i];
			if (energy > max_value) {
				max_value = energy;
			}
			if (energy < min_value) {
				min_value = energy;
			}
		}

		var	shift = 0.0 - min_value;
		var	range = max_value - min_value;

		for ( var i = 0; i < this.coords.length; ++i) {
			var energy = energy_vectors[0][i];
			var hue = (energy + shift) / range;
			this.nodes.colors[i].setHSV(hue, 1.0, 1.0); //*
		}
	}
	this.add_edges = function() {
		var nodes = this.nodes;
		var edges = this.edges;
		var lines = this.lines;

		var geometry = new THREE.Geometry();
		var colours = [];

		for ( var i = 0; i < edges.length; ++i) {
			var a = new THREE.Vertex( new THREE.Vector3( 0.0, 0.0, 0.0) )
			var b = new THREE.Vertex( new THREE.Vector3( 0.0, 0.0, 0.0) )
			geometry.vertices.push( a );
			geometry.vertices.push( b );
		}

		for (var i =0; i < geometry.vertices.length; ++i) {
			colours[ i ] = new THREE.Color( 0xffffff );
			colours[ i ].setHSV( i / geometry.vertices.length, 1.0, 1.0 );
		}

		geometry.colors = colours;

		var line_material = new THREE.LineBasicMaterial({
			color: 0xffffff,
			opacity: 0.5,
			linewidth: 2
		} );
		line_material.vertexColors = true;

		var line = new THREE.Line( geometry, line_material, THREE.LinePieces );
		line.scale.x = line.scale.y = line.scale.z = 1.0;

		line.updateMatrix();
		line.geometry.__dirtyVertices = true;
		line.geometry.__dirtyColors = true;
		this.add_to_scene([line]);
		this.lines = line;
		this.move_edges();

	}
	this.update_edge_colours = function() {
		var edges = this.edges;
		var energies = this.energies[0];
		var colours = this.lines.geometry.colors;

		var min_value = 1e9;
		var max_value = -1e9;


		for ( var i = 0; i < edges.length -1; ++i) {
			var u = edges[i][0];
			var v = edges[i][1];

			var e_u = energies[u];
			var e_v = energies[v];
			var diff = Math.abs(e_u - e_v);
			var energy = diff;
			if (energy > max_value) {
				max_value = energy;
			}
			if (energy < min_value) {
				min_value = energy;
			}
		}

		var	shift = 0.0 - min_value;
		var	range = max_value - min_value;

		for ( var i = 0; i < edges.length -1; ++i) {
			var u = edges[i][0];
			var v = edges[i][1];

			var e_u = energies[u];
			var e_v = energies[v];
			var diff = Math.abs(e_u - e_v);
			var hue = (diff + shift) / range;
			console.log(hue);
			colours[2*i].setRGB(hue, hue, hue);
			colours[2*i+1].setRGB(hue, hue, hue);
		}
		this.lines.geometry.__dirtyColors = true;
	}
	this.move_edges = function() {
		var nodes = this.nodes;
		var edges = this.edges;
		var lines = this.lines;

		for ( var i = 0; i < edges.length -1; ++i) {
			var u = edges[i][0];
			var v = edges[i][1];
			var a = nodes.vertices[u];
			var b = nodes.vertices[v];
			lines.geometry.vertices[2*i].position.x = a.position.x;
			lines.geometry.vertices[2*i].position.y = a.position.y;
			lines.geometry.vertices[2*i].position.z = a.position.z;

			lines.geometry.vertices[2*i+1].position.x = b.position.x;
			lines.geometry.vertices[2*i+1].position.y = b.position.y;
			lines.geometry.vertices[2*i+1].position.z = b.position.z;
		}
		lines.geometry.__dirtyVertices = true;

	}
	this.add_to_scene = function(objects) {
		var scope = this;
		$.each(objects, function(index, value) {
			scope.scene.addObject( value );
		});
	}
}