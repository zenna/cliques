function Graph(data, scene) {
	this.coords = data.coords;
	this.energies = data.energies;
	this.edges = data.edges;
	this.mean_position = {
		x:0.0,
		y:0.0,
		z:0.0
	};
	this.scene = scene;
	this.materials = [], this.nodes = [], this.lines = [];
	this.num_dim = this.coords[0].length;

	this.place_nodes = function(dim1, dim2, dim3) {
		var materials = this.materials;
		var nodes = this.nodes;
		var coords = this.coords;
		var mean_position = this.mean_position;

		// for ( var i = 0; i < coords.length; ++i) {
		// var vector = new THREE.Vector3(0.0,0.0,0.0);
		// var geometry = new THREE.Geometry();
		// geometry.vertices.push(new THREE.Vertex(vector));
		//
		// var hue  = 0.1 * (i % 9);
		// var color = [hue, 1.0, 1.0];
		// var material = new THREE.ParticleBasicMaterial({
		// size:40
		// });
		// material.color.setHSV(color[0], color[1], color[2]);
		// var node = new THREE.ParticleSystem(geometry, material);
		// nodes.push(node);
		// materials.push(material);
		// }

		var geometry = new THREE.Geometry();
		var colours = [];

		for ( var i = 0; i < coords.length; ++i) {
			x = 2000 * Math.random() - 1000;
			y = 2000 * Math.random() - 1000;
			z = 2000 * Math.random() - 1000;
			vector = new THREE.Vector3( x, y, z );
			geometry.vertices.push( new THREE.Vertex( vector ) );
			colours[i] = new THREE.Color( 0xffffff );
			var hue  = 0.1 * (i % 9);
			colours[ i ].setHSV(hue, 1.0, 1.0 );
		}

		// sprite = THREE.ImageUtils.loadTexture( "textures/sprites/ball.png" );//
		geometry.colors = colours;

		var material = new THREE.ParticleBasicMaterial({
			size: 40,
			vertexColors: true
		} );
		material.color.setHSV( 1.0, 0.2, 0.8 );

		var particles = new THREE.ParticleSystem( geometry, material );
		particles.sortParticles = true;
		particles.updateMatrix();
		this.nodes = geometry;
		this.particles = particles;

		this.move_nodes();
		this.add_to_scene([particles]);
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
		this.move_edges();
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
	this.add_edges = function(edges_vectors) {
		var nodes = this.nodes;
		var edges = this.edges;
		var lines = this.lines;

		line_material = new THREE.LineBasicMaterial({
			color: 0xffaa00,
			opacity: 0.5,
			linewidth: 2
		} );

		var line = new THREE.Line( this.nodes, line_material, THREE.LinePieces );
		var colour = new THREE.Color( 0xffffff );
		line.colors = colour;
		line.updateMatrix();
		this.line = line;
		line.geometry.__dirtyVertices = true;

// 
		// for (var i=0;i< edges.length; ++i) {
			// var u = edges[i][0];
			// var v = edges[i][1];
			// console.log(u + "-" + v);
			// var line_geometry = new THREE.Geometry();
			// line_geometry.vertices.push(nodes[u].geometry.vertices[0]);
			// line_geometry.vertices.push(nodes[v].geometry.vertices[0]);
			// colors = [];
// 
			// var position, index;
			// for (var j = 0; j < line_geometry.vertices.length; j ++ ) {
// 
				// colors[ j ] = new THREE.Color( 0xffffff );
				// colors[ j ].setHSV( 0.6, 1.0, 1.0 );
			// }
			// line_geometry.colors = colors;
// 
			// line_material = new THREE.LineBasicMaterial({
				// color: 0xffffff,
				// opacity: 0.5,
				// linewidth: 2
			// } );
			// line_material.vertexColors = true;
			// var line = new THREE.Line( line_geometry, line_material );
			// lines.push(line);
		// }

		this.add_to_scene([line]);

	}
	this.move_edges = function() {
		var nodes = this.nodes;
		var edges = this.edges;
		var lines = this.lines;

		for ( var i = 0; i < lines.length; ++i) {
			var u = edges[i][0];
			var v = edges[i][1];
			lines[i].geometry.vertices[0].position.x = nodes[u].geometry.vertices[0].position.x;
			lines[i].geometry.vertices[0].position.y = nodes[u].geometry.vertices[0].position.y;
			lines[i].geometry.vertices[0].position.z = nodes[u].geometry.vertices[0].position.z;

			lines[i].geometry.vertices[1].position.x = nodes[v].geometry.vertices[0].position.x;
			lines[i].geometry.vertices[1].position.y = nodes[v].geometry.vertices[0].position.y;
			lines[i].geometry.vertices[1].position.z = nodes[v].geometry.vertices[0].position.z;
			lines[i].geometry.__dirtyVertices = true;
		}
	}
	this.add_to_scene = function(objects) {
		var scope = this;
		$.each(objects, function(index, value) {
			scope.scene.addObject( value );
		});
	}
}