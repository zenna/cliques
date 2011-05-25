function setup() {
	// set the scene size
	var WIDTH = 400, HEIGHT = 300;

	// set some camera attributes
	var VIEW_ANGLE = 45, ASPECT = WIDTH / HEIGHT, NEAR = 0.1, FAR = 10000;

	// get the DOM element to attach to
	// - assume we've got jQuery to hand
	container = $('#container');

	// create a WebGL renderer, camera
	// and a scene
	var renderer = new THREE.WebGLRenderer();
	var camera = new THREE.Camera(VIEW_ANGLE, ASPECT, NEAR, FAR);
	var scene = new THREE.Scene();

	// the camera starts at 0,0,0 so pull it back
	camera.position.z = 300;

	// start the renderer
	renderer.setSize(WIDTH, HEIGHT);

	// attach the render-supplied DOM element
	container.append(renderer.domElement);

	// create the sphere's material
	var sphereMaterial = new THREE.MeshLambertMaterial({
		color : 0xCC0000
	});
	// set up the sphere vars
	var radius = 50, segments = 16, rings = 16;

	// create a new mesh with sphere geometry -
	// we will cover the sphereMaterial next!
	sphere = new THREE.Mesh(new THREE.Sphere(radius, segments, rings),
	sphereMaterial);

	sphere.scale.x = 2;

	// add the sphere to the scene
	scene.addChild(sphere);

	// create a point light
	var pointLight = new THREE.PointLight(0xFFFFFF);

	// set its position
	pointLight.position.x = 10;
	pointLight.position.y = 50;
	pointLight.position.z = 130;

	// add to the scene
	scene.addLight(pointLight);

	// draw!
	renderer.render(scene, camera);
}

function setup2() {

	if (!Detector.webgl)
		Detector.addGetWebGLMessage();

	var container, stats;
	var camera, scene, renderer, particles, geometry, materials = [], parameters, i, h, color;
	var mouseX = 0, mouseY = 0;

	var windowHalfX = window.innerWidth / 2;
	var windowHalfY = window.innerHeight / 2;

	init();
	animate();

	function init() {

		container = document.createElement('div');
		document.body.appendChild(container);

		camera = new THREE.Camera(75, window.innerWidth / window.innerHeight,
		1, 3000);
		camera.position.z = 1000;

		scene = new THREE.Scene();
		scene.fog = new THREE.FogExp2(0x000000, 0.0007);

		geometry = new THREE.Geometry();

		for ( var i = 0; i < vectors.length; ++i) {
			vector = new THREE.Vector3(vectors[i][0] * 2000 - 1000,
			vectors[i][1] * 2000 - 1000, vectors[i][2] * 2000 - 1000);
			geometry.vertices.push(new THREE.Vertex(vector));
		}

		parameters = [ [ [ 1.0, 1.0, 1.0 ], 5 ], [ [ 0.95, 1, 1 ], 4 ],
		[ [ 0.90, 1, 1 ], 3 ], [ [ 0.85, 1, 1 ], 2 ],
		[ [ 0.80, 1, 1 ], 1 ] ];
		// parameters = [ [ 0xff0000, 5 ], [ 0xff3300, 4 ], [ 0xff6600, 3 ], [
		// 0xff9900, 2 ], [ 0xffaa00, 1 ] ];
		// parameters = [ [ 0xffffff, 5 ], [ 0xdddddd, 4 ], [ 0xaaaaaa, 3 ], [
		// 0x999999, 2 ], [ 0x777777, 1 ] ];

		for (i = 0; i < parameters.length; i++) {

			size = parameters[i][1];
			color = parameters[i][0];

			// materials[i] = new THREE.ParticleBasicMaterial( { color: color,
			// size: size } );

			materials[i] = new THREE.ParticleBasicMaterial({
				size : size
			});
			materials[i].color.setHSV(color[0], color[1], color[2]);

			particles = new THREE.ParticleSystem(geometry, materials[i]);
			// particles.rotation.x = Math.random() * 6;
			// particles.rotation.y = Math.random() * 6;
			// particles.rotation.z = Math.random() * 6;
			scene.addObject(particles);

		}

		renderer = new THREE.WebGLRenderer({
			antialias : false
		});
		renderer.setSize(window.innerWidth, window.innerHeight);
		container.appendChild(renderer.domElement);

		stats = new Stats();
		stats.domElement.style.position = 'absolute';
		stats.domElement.style.top = '0px';
		container.appendChild(stats.domElement);

		document.addEventListener('mousemove', onDocumentMouseMove, false);
		document.addEventListener('touchstart', onDocumentTouchStart, false);
		document.addEventListener('touchmove', onDocumentTouchMove, false);

	}

	function onDocumentMouseMove(event) {
		mouseX = event.clientX - windowHalfX;
		mouseY = event.clientY - windowHalfY;
	}

	function onDocumentTouchStart(event) {
		if (event.touches.length == 1) {
			event.preventDefault();
			mouseX = event.touches[0].pageX - windowHalfX;
			mouseY = event.touches[0].pageY - windowHalfY;
		}
	}

	function onDocumentTouchMove(event) {
		if (event.touches.length == 1) {
			event.preventDefault();
			mouseX = event.touches[0].pageX - windowHalfX;
			mouseY = event.touches[0].pageY - windowHalfY;
		}
	}

	function animate() {
		requestAnimationFrame(animate);
		render();
		stats.update();
	}

	function render() {
		var time = new Date().getTime() * 0.00005;
		camera.position.x += (mouseX - camera.position.x) * 0.05;
		camera.position.y += (-mouseY - camera.position.y) * 0.05;
		for (i = 0; i < scene.objects.length; i++) {
			scene.objects[i].rotation.y = time * (i < 4 ? i + 1 : -(i + 1));
		}

		for (i = 0; i < materials.length; i++) {
			color = parameters[i][0];
			h = (360 * (color[0] + time) % 360) / 360;
			materials[i].color.setHSV(h, color[1], color[2]);

		}
		renderer.render(scene, camera);
	}
}