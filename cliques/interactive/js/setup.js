app = {
	points:[],

	setup: function(vectors) {
		ok = 10;
		this.vectors = vectors;

		if (!Detector.webgl)
			Detector.addGetWebGLMessage();

		var container, stats;
		var camera, scene, renderer, particles, geometry, parameters, i, h, color;
		var mouseX = 0, mouseY = 0;
		materials = [];

		var windowHalfX = window.innerWidth / 2;
		var windowHalfY = window.innerHeight / 2;
		parent = this;

		this.init();
		animate();

		function animate() {
			requestAnimationFrame( animate );
			parent.render();
			parent.stats.update();
		}

	},
	init: function() {
		var vectors = this.vectors;
		var mean_position = {
			x:0.0,
			y:0.0,
			z:0.0
		};

		container = document.createElement('div');
		document.body.appendChild(container);

		scene = new THREE.Scene();
		this.scene = scene;

		for ( var i = 0; i < vectors.length; ++i) {
			vector = new THREE.Vector3(vectors[i][0] * 2000 - 1000,
			vectors[i][1] * 2000 - 1000, vectors[i][2] * 2000 - 1000);
			var geometry = new THREE.Geometry();
			geometry.vertices.push(new THREE.Vertex(vector));

			var hue  = 0.1 * (i % 9);
			var color = [hue, 1.0, 1.0];
			material = new THREE.ParticleBasicMaterial({
				size:40
			});
			material.color.setHSV(color[0], color[1], color[2]);
			particles = new THREE.ParticleSystem(geometry, material);
			scene.addObject(particles);
			materials.push(material);

			mean_position.x = mean_position.x + vector.x;
			mean_position.y = mean_position.y + vector.y;
			mean_position.z = mean_position.z + vector.z;
			this.points.push(geometry);
		}

		mean_position.x = mean_position.x / vectors.length;
		mean_position.y = mean_position.y / vectors.length;
		mean_position.z = mean_position.z / vectors.length;

		this.renderer = new THREE.WebGLRenderer({
			antialias : true
		});
		var renderer = this.renderer;
		renderer.setSize(window.innerWidth, window.innerHeight);
		container.appendChild(renderer.domElement);

		this.stats = new Stats();
		stats = this.stats;
		stats.domElement.style.position = 'absolute';
		stats.domElement.style.top = '0px';
		container.appendChild(stats.domElement);

		camera = new THREE.TrackballCamera({

			fov: 75,
			aspect: window.innerWidth / window.innerHeight,
			near: 1,
			far: 10000,

			rotateSpeed: 1.0,
			zoomSpeed: 1.2,
			panSpeed: 0.2,

			noZoom: false,
			noPan: false,

			staticMoving: false,
			dynamicDampingFactor: 0.3,

			keys: [ 65, 83, 68 ], // [ rotateKey, zoomKey, panKey ],

			domElement: renderer.domElement,

		});

		camera.position.z = 4000;
		camera.target.position.x = mean_position.x;
		camera.target.position.y = mean_position.y;
		camera.target.position.z = mean_position.z;
	},
	render: function() {
		for (var i=0; i<scene.objects.length;++i) {
			scene.objects[i].geometry.vertices[0].position.x += 1.0;
			scene.objects[i].geometry.vertices[0].position.y += 1.0;
			scene.objects[i].geometry.vertices[0].position.z += 1.0;
			scene.objects[i].geometry.__dirtyVertices = true;
		}
		$(document).trigger('render');
		this.renderer.render(scene, camera);
	},
	add_colours: function(energy_vectors) {
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

		for ( var i = 0; i < this.vectors.length; ++i) {
			var energy = energy_vectors[0][i];
			var hue = (energy + shift) / range;
			materials[i].color.setHSV(hue, 1.0, 1.0);
		}
	},
	add_edges: function(edges_vectors) {
		var points = this.points;
		var scene = this.scene;

		for (var i=0;i< edges_vectors.length; ++i) {
			var u = edges_vectors[i][0];
			var v = edges_vectors[i][1];
			console.log(u + "-" + v);
			var line_geometry = new THREE.Geometry();
			line_geometry.vertices.push(points[u].vertices[0]);
			line_geometry.vertices.push(points[v].vertices[0]);
			colors = [];

			var position, index;
			for (var j = 0; j < line_geometry.vertices.length; j ++ ) {

				colors[ j ] = new THREE.Color( 0xffffff );
				colors[ j ].setHSV( 0.6, 1.0, 1.0 );
			}
			line_geometry.colors = colors;

			line_material = new THREE.LineBasicMaterial({
				color: 0xffffff,
				opacity: 0.5,
				linewidth: 2
			} );
			line_material.vertexColors = true;
			var line = new THREE.Line( line_geometry, line_material );
			this.scene.addObject( line );
		}
	}
}