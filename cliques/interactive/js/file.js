$(document).ready(function() {
    function handleFileSelect(evt) {
        vectors = []
        var files = evt.target.files; // FileList object
        f = files[0];
        if (f) {
            var reader = new FileReader();
            reader.onload = function(e) {
                coords = e.target.result.split("\n");
                for (var i=0;i<coords.length;++i) {
                    var vector = coords[i].split(" ");
                    var float_vector = [];
                    for (var j=0;j<vector.length;++j) {
                        if (vector[j] != " " && vector[j] != "") {
                            float_vector.push(parseFloat(vector[j]));
                        } 
                    }
                    vectors.push(float_vector);
                }
                console.log(e.target.result);
                setup2();
            }
            reader.readAsText(f);
        }
        else {
            alert("no file");
        }
    }
    
    function handle_energy_files(evt) {
        energy_vectors = []
        var files = evt.target.files; // FileList object
        f = files[0];
        if (f) {
            var reader = new FileReader();
            reader.onload = function(e) {
                coords = e.target.result.split("\n");
                for (var i=0;i<coords.length;++i) {
                    var vector = coords[i].split(" ");
                    var float_vector = [];
                    for (var j=0;j<vector.length;++j) {
                        if (vector[j] != " " && vector[j] != "") {
                            float_vector.push(parseFloat(vector[j]));
                        } 
                    }
                    energy_vectors.push(float_vector);
                }
                add_colours();
                //$(document).bind('render', { 'bar' : 'bam' }, add_colour);
            }
            reader.readAsText(f);
        }
        else {
            alert("no file");
        }
    }

    document.getElementById('files').addEventListener('change',
        handleFileSelect, false);
        
    document.getElementById('energy_files').addEventListener('change',
        handle_energy_files, false);
    })