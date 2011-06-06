function setup2() {

	if (!Detector.webgl)
		Detector.addGetWebGLMessage();

	var container, stats;
	var camera, scene, renderer, particles, geometry, parameters, i, h, color;
	var mouseX = 0, mouseY = 0;
	materials = [];

	var windowHalfX = window.innerWidth / 2;
	var windowHalfY = window.innerHeight / 2;

	init();
	animate();

	function init() {

		var mean_position = {
			x:0.0,
			y:0.0,
			z:0.0
		};

		container = document.createElement('div');
		document.body.appendChild(container);

		scene = new THREE.Scene();
		//scene.fog = new THREE.FogExp2(0x000000, 0.0007);

		// geometry = new THREE.Geometry();	

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
		}

		mean_position.x = mean_position.x / vectors.length;
		mean_position.y = mean_position.y / vectors.length;
		mean_position.z = mean_position.z / vectors.length;

		// var color = [1.0, 1.0, 1.0];
		// material = new THREE.ParticleBasicMaterial({
			// size:40
		// });
		// material.color.setHSV(color[0], color[1], color[2]);
// 
		// for ( var i = 0; i < vectors.length; ++i) {
			// var hue  = 0.1 * (i % 9);
			// var color = [hue, 1.0, 1.0];
			// material = new THREE.ParticleBasicMaterial({
				// size:40
			// });
			// material.color.setHSV(color[0], color[1], color[2]);
			// materials.push(material);
		// }
// 
		// particles = new THREE.ParticleSystem(geometry, materials);
		// scene.addObject(particles);

		renderer = new THREE.WebGLRenderer({
			antialias : false
		});
		renderer.setSize(window.innerWidth, window.innerHeight);
		container.appendChild(renderer.domElement);

		stats = new Stats();
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
	}

	function animate() {
		requestAnimationFrame(animate);
		render();
		stats.update();
	}

	function render() {
		var time = new Date().getTime() * 0.00005;
		/*for (i = 0; i < scene.objects.length; i++) {
		scene.objects[i].rotation.y = time * (i < 4 ? i + 1 : -(i + 1));
		}*/

		// for (i = 0; i < materials.length; i++) {
		// color = parameters[i][0];
		// h = (360 * (color[0] + time) % 360) / 360;
		// materials[i].color.setHSV(h, color[1], color[2]);
		//
		// }
		$(document).trigger('render');
		renderer.render(scene, camera);
	}
}

function add_colours() {
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
	
	for ( var i = 0; i < vectors.length; ++i) {
		var energy = energy_vectors[0][i];
		var hue = (energy + shift) / range;
		materials[i].color.setHSV(hue, 1.0, 1.0);
	}
}

function add_colour() {
	for (i = 0; i < energy_vectors.length; i++) {
		console.log(i);
		h = (360 * (color[0] + time) % 360) / 360;
		material.color.setHSV(color[0], color[1], color[2]);
	}
}