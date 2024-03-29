function App(width, height, name) {
	var	width = typeof width  == "undefined" ? window.innerWidth : width;
	var	height = typeof height  == "undefined" ? window.innerHeight : height;
	var name = typeof name  == "undefined" ? "" : name;

	this.setup = function() {
		if (!Detector.webgl)
			Detector.addGetWebGLMessage();

		var parent = this;

		this.init(name);
		animate();

		function animate() {
			requestAnimationFrame( animate );
			parent.render();
		}

	}
	this.init = function(name) {
		var container = document.createElement('div');
		document.body.appendChild(container);
		this.scene = new THREE.Scene();
		this.renderer = new THREE.WebGLRenderer({
			antialias : true
		});
		var renderer = this.renderer;
		renderer.setSize(width, height);
		container.appendChild(renderer.domElement);
		$(renderer.domElement).addClass(name);

		this.camera = new THREE.TrackballCamera({

			fov: 75,
			aspect: window.innerWidth / window.innerHeight,
			near: 1,
			far: 50000,

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
		document.addEventListener( 'click', this.onDocumentMouseClick, false );
		document.addEventListener( 'mousemove', this.hover, false );

	}
	this.hover = function(event) {
		var context = landscape_view.renderer.context;
		var mouseX = event.clientX;
		var mouseY = event.clientY;
		var invertedMouseY = $(landscape_view.renderer.domElement).height() - mouseY;alert

		landscape_view.renderer.render(landscape_view.scene, landscape_view.camera);
		var arr = new Uint8Array(4);
		context.readPixels(mouseX, invertedMouseY, 1, 1, context.RGBA, context.UNSIGNED_BYTE, arr);
		var id = cliques.rgbToInt(arr);
		if (id != 0) {
			document.body.style.cursor = 'pointer';
		} else {
			document.body.style.cursor = 'default';
		}
		$('.processPopup').css('left',mouseX+15).css('top',mouseY+15);


	}
	this.onDocumentMouseClick = function(event) {
		var app = landscape_view;
		var renderer = app.renderer;
		var render_target = new THREE.WebGLRenderTarget( window.innerWidth, window.innerHeight, {
			minFilter: THREE.LinearFilter,
			magFilter: THREE.NearestFilter,
			format: THREE.RGBFormat
		});

		var original_colours = landscape.nodes.colors;
		landscape.nodes.colors = landscape.node_id_colours;
		var opacity = landscape.lines.materials[0].opacity;
		landscape.lines.materials[0].opacity = 0.0;
		landscape.nodes.__dirtyColors = true;
		// var id_scene = new Three.scene()
		// id_scene.addObject( lanscape.nodes );

		renderer.render( app.scene, app.camera, render_target, true );
		// Return to normal
		landscape.nodes.colors = original_colours;
		landscape.lines.materials[0].opacity = opacity;

		var mouseX = event.clientX;
		var mouseY = event.clientY;
		// readPixels returns coords with bottom_left as origin
		// whereas clientX/Y have top_left origin
		var invertedMouseY = $(app.renderer.domElement).height() - mouseY;

		var context = renderer.context;
		var arr = new Uint8Array(4);
		context.readPixels(mouseX, invertedMouseY, 1, 1, context.RGBA, context.UNSIGNED_BYTE, arr);
		var id = cliques.rgbToInt(arr)
		console.log(arr[0], arr[1],arr[2],id);

		var partition = landscape.partitions[id];
		if (typeof partition == 'undefined') {
			return;
		}
		console.log(landscape.partitions[id]);
		var colorMap;
		if (landscape.landscapeType == 'community') {
			colorMap = new cliques.CommunityColorMap()
		} else {
			colorMap = new cliques.PartitionColorMap()
		}

		var partition_colors = colorMap.getColors(partition);
		orig_graph.paint_nodes(partition_colors);
		orig_graph.match_edge_colours_to_node();
		
		//Quick way is just to look up stability in process toolbox and then find stability and render an insibile box
		// More elegant way is to let each process influence a box shown by the cusor with metadata. e.g. if we also want to
		// show metadata for a node, e.g. its probability of going to all the different basins.
		// How to do this, when we click a node and have ound its id, we want to trigger off all the basin.  
		// Complexation
		toolbox.popupNodeDataFunctions(id);
		
	}
	this.render = function() {
		$(document).trigger('render');
		this.renderer.render(this.scene, this.camera);

	}
	this.save = function(sizeX, sizeY) {
		this.renderer.setSize(sizeX, sizeY);
		this.render();
		window.open(this.renderer.domElement.toDataURL("image/png"));
	}
}