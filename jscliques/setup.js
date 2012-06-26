$(document).ready(function() {
    // var jqxkr = $.getJSON('js/data/tree.json', function(data) {
    // load_landscape(data);
    // }).error( function(data) {
    // alert("error, could not load");
    // })
    //
    function load_landscape(data) {
        var width = parseInt($("#xLength").val());
        var height = parseInt($("#yLength").val());
        landscape_view = new App(width, height, "landscapeView");
        landscape_view.setup();

        // Create a new graph for every graph in data_set
        landscapes = [];
        for (var i=0;i<data.length;++i) {
            var landscape_layer = new Graph(data[i], landscape_view.scene, 'landscape');
            landscape_layer.create_random_offsets(0.02);
            landscape_layer.place_nodes();
            landscape_layer.add_edges({
                opacity : 0.6
            });
            landscapes.push(landscape_layer);
        }
        landscape = landscapes[0];
        landmas = 21;
        // landscape = new Graph(data, landscape_view.scene, 'landscape');
        // landscape.create_random_offsets(0.02);
        // landscape.place_nodes();

        // Position tthe camera to target the centre of the 0th landscape
        landscape_view.camera.target.position.x = landscapes[0].mean_position.x;
        landscape_view.camera.target.position.y = landscapes[0].mean_position.y;
        landscape_view.camera.target.position.z = landscapes[0].mean_position.z;
        landscape_view.camera.position.z = 4000;

        // landscape.add_edges({
        //     opacity : 0.6
        // });
        //landscape.update_edge_colours();

        graph_view = new App(300, 300, 'graph_view');
        graph_view.setup();
        orig_graph = new Graph(data[0].graph, graph_view.scene);
        orig_graph.create_random_offsets(0);
        orig_graph.place_nodes();
        graph_view.camera.target.position.x = orig_graph.mean_position.x;
        graph_view.camera.target.position.y = orig_graph.mean_position.y;
        graph_view.camera.target.position.z = orig_graph.mean_position.z;
        graph_view.camera.position.z = 1300;
        orig_graph.add_edges({
            opacity : 1.0
        });
        box = new cliques.SelectionBox(graph_view.camera.domElement);
        toolbox = new ProcessToolbox(landscape_view);

        // Add Processes - HACK, just do for 0th landscape
        for(var i = 0; i < data[0].processes.length; ++i) {
            var process = data[0].processes[i];
            switch(process.type) {
                case 'basins':
                    toolbox.addProcess(new BasinProcess(process, landscapes[0]));
                    break;
                case 'stability':
                    toolbox.addProcess(new NodeProcess(process, landscapes[0]));
                    break;
            }
        }
    }

    function saveImage() {
        var sizeX = parseInt($("#xLength").val());
        var sizeY = parseInt($("#yLength").val());
        landscape_view.save(sizeX, sizeY);
        //imageElement.src = myImage;
        // return dataURL.replace(/^data:image\/(png|jpg);base64,/, "");
        // https://developer.mozilla.org/En/HTML/Canvas/Pixel_manipulation_with_canvas
    }

    // For handling button to load file
    function handleFileSelect(evt) {
        var files = evt.target.files;
        // FileList object
        f = files[0];
        if(f) {
            var reader = new FileReader();
            reader.onload = function(e) {
                var landscape_json = $.parseJSON(e.target.result);
                load_landscape(landscape_json);
            }


            reader.readAsText(f);
        }
        else {
            alert("no file");
        }
    }


    document.getElementById('files').addEventListener('change', handleFileSelect, false);

    $('#save').click(function() {
        saveImage();
    })


    $('#genUrl').click(function() {
        var camera_data = {};
        camera_data['px'] = landscape_view.camera.position.x;
        camera_data['py'] = landscape_view.camera.position.y;
        camera_data['pz'] = landscape_view.camera.position.z;
        camera_data['tx'] = landscape_view.camera.target.position.x;
        camera_data['ty'] = landscape_view.camera.target.position.y;
        camera_data['tz'] = landscape_view.camera.target.position.z;
        var cameraDataString = $.param(camera_data);
        alert(cameraDataString);
    })


    $(window).hashchange(function() {
        var hash = location.hash;
        landscape_view.camera.position.x = parseFloat($.url(hash).fparam('px'));
        landscape_view.camera.position.y = parseFloat($.url(hash).fparam('py'));
        landscape_view.camera.position.z = parseFloat($.url(hash).fparam('pz'));
        landscape_view.camera.target.position.x = parseFloat($.url(hash).fparam('tx'));
        landscape_view.camera.target.position.y = parseFloat($.url(hash).fparam('ty'));
        landscape_view.camera.target.position.z = parseFloat($.url(hash).fparam('tz'));
    })

    function loadHashData(hash) {

    }

    var screenWidth = typeof width == "undefined" ? window.innerWidth : width;
    var screenHeight = typeof height == "undefined" ? window.innerHeight : height;
    $("#xLength").val(screenWidth);
    $("#yLength").val(screenHeight);

    // On hash change, modify program

})