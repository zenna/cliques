app = {
	setup: function() {
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
		container = document.createElement('div');
		document.body.appendChild(container);
		this.scene = new THREE.Scene();
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

		this.camera = new THREE.TrackballCamera({

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
// 
		// camera.position.z = 4000;
		// camera.target.position.x = mean_position.x;
		// camera.target.position.y = mean_position.y;
		// camera.target.position.z = mean_position.z;
	},
	render: function() {
		// for (var i=0; i<this.scene.objects.length;++i) {
			// scene.objects[i].geometry.vertices[0].position.x += 1.0;
			// scene.objects[i].geometry.vertices[0].position.y += 1.0;
			// scene.objects[i].geometry.vertices[0].position.z += 1.0;
			// scene.objects[i].geometry.__dirtyVertices = true;
		// }
		$(document).trigger('render');
		this.renderer.render(this.scene, this.camera);
	},
}