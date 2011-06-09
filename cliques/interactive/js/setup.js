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

			//domElement: renderer.domElement,

		});

		document.addEventListener( 'click', this.onDocumentMouseClick, false );
		//
		// camera.position.z = 4000;
		// camera.target.position.x = mean_position.x;
		// camera.target.position.y = mean_position.y;
		// camera.target.position.z = mean_position.z;
	},
	onDocumentMouseClick: function(event) {
		var renderer = app.renderer;
		render_target = new THREE.WebGLRenderTarget( window.innerWidth, window.innerHeight, {
			minFilter: THREE.LinearFilter,
			magFilter: THREE.NearestFilter,
			format: THREE.RGBFormat
		});

		var original_colours = landscape.nodes.colors;
		landscape.nodes.colors = landscape.node_id_colours;
		renderer.render( app.scene, app.camera, render_target, true );
		landscape.nodes.colors = original_colours;
		
		var mouseX = event.clientX;
		var mouseY = event.clientY;
		// readPixels returns coords with bottom_left as origin
		// whereas clientX/Y have top_left origin
		var invertedMouseY = $(app.renderer.domElement).height() - mouseY;

		var context = renderer.context;
		var arr = new Uint8Array(4);
		context.readPixels(mouseX, invertedMouseY, 1, 1, context.RGBA, context.UNSIGNED_BYTE, arr);
		console.log(arr[0],arr[1],arr[2]);
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