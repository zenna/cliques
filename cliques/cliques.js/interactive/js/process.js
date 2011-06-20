louvain = {
	type:'louvain',
	dependencies: {
		type:'stability',
		time:0.234
	},
	data:[{
		partition:'034202',
		action:'consider',
		time:0
	},{
		partition:'034202',
		action:'consider',
		time:0
	},{
		partition:'034202',
		action:'position',
		time:1
	},
	]
}

louvain2 = {
	type:'louvain',
	dependencies: {
		type:'stability',
		time:4.3
	},
	data:[{
		partition:'034202',
		action:'consider',
		time:0
	},{
		partition:'034202',
		action:'consider',
		time:0
	},{
		partition:'034202',
		action:'position',
		time:1
	},
	]
}

var NodeProcess = function(process, graph, dynamicColors) {
	this.type = process.type;
	this.data = process.data;
	this.graph = graph;

	// if true, scale coloring to data at particular point
	// as opposed to absolute range over all data points
	var dynamicColors = true;
	this.dynamicColors = dynamicColors;

	// Initialse range of scalers in color maps to data
	if (dynamicColors == true) {
		this.colorMaps = [];
		for (var i=0;i<this.data.length;++i) {
			var colorMap = new cliques.EnergyColorMap();
			for (var j = 0; j < this.data[i].values.length;++j) {
				colorMap.updateExtrema(this.data[i].values[j]);
			}
			this.colorMaps.push(colorMap);
		}
	}

	this.colorMap = new cliques.EnergyColorMap();
	for (var i=0;i<this.data.length;++i) {
		for (var j = 0; j < this.data[i].values.length;++j) {
			this.colorMap.updateExtrema(this.data[i].values[j]);
		}
	}
};
NodeProcess.prototype.render = function(dataId) {
	console.log(dataId);
	var rgbs = [];
	var data = this.data[dataId];
	if (this.dynamicColors == true) {
		for (var i=0;i<data.values.length;++i) {
			var color = this.colorMaps[dataId].getColor(data.values[i]);
			rgbs.push(color);
		}
	} else {
		for (var i=0;i<data.values.length;++i) {
			var color = this.colorMap.getColor(data.values[i]);
			rgbs.push(color);
		}
	}
	this.graph.paint_nodes(rgbs);
	this.graph.match_edge_colours_to_node();
	var norm_values = [];
	for (var i=0;i<data.values.length;++i) {
		norm_values.push(this.colorMaps[dataId].scaler.scaleValue(data.values[i]) * 2.0);
	}
	this.graph.move_nodes(norm_values, 0, 1);
}
NodeProcess.prototype.hide = function(time) {

}
var BasinProcess = function(process, graph) {
	this.type = process.type;
	this.data = process.data;
	this.graph = graph;
}
BasinProcess.prototype.render = function(dataId) {
	var values = this.data[dataId].values;
	var colorMap = new cliques.EnergyColorMap();
	colorMap.updateExtrema(0);
	colorMap.updateExtrema(values.length);
	var rgbs = [];

	var norm_values = [];
	for (var i =0;i<values.length;++i) {
		var basinColor = colorMap.getColor(i);
		// graph.highlight_node(values['basin'], basinColor));
		for (var j = 0;j<values[i].nodes.length;++j) {
			var node_id = values[i].nodes[j];
			rgbs[node_id] = basinColor;
			norm_values[node_id] = colorMap.scaler.scaleValue(i);
		}
		rgbs[values[i]['basin']] = colorMap.getColor((i+1)*0.6);

	}
	
	this.graph.paint_nodes(rgbs);
	// this.graph.move_nodes(norm_values, 0, 1);
	// this.graph.match_edge_colours_to_node();


	// var alpha = {
	// type:'basin',
	// dependencies:['stability'],
	// data:[{time:0.234, values:[{basin:0,nodes:[0,2,3]},
	// {basin:1,nodes:[3,4,5]}];
	// this.graph.match_edge_colours_to_node();
	// }

}
var ProcessToolbox = function(view) {
	this.domElement = $("<div></div>").addClass('toolbox');
	$(view.renderer.domElement).parent().append(this.domElement);
	this.processes = {};
}
ProcessToolbox.prototype.addProcess = function(process) {
	var type = process.type;
	if (typeof this.processes[type] == "undefined") {
		this.processes[type] = [];
	}
	this.processes[type].push(process);
	this.updateDomElement();
}
dynamicProcessTemplate = "<div><span class='processName'></span>@\
<input class = 'timeSlider' name='r' type='range' min='0' max='311' value='0'>\
<span class='value'></span>\
</div>\
<br/>\
"

processSliderTemplate = "<input class = 'processSlider' name='r' type='range' min='0' max='10' value='0'>"

ProcessToolbox.prototype.updateDomElement = function() {
	// Generate html from processes

	this.domElement.html('');
	for (var process in this.processes) {
		var processGroup = this.processes[process];
		var processHtml = $(dynamicProcessTemplate);
		var timeSlider = processHtml.children('input.timeSlider')
		this.updateSliderFromProcess(processGroup[0], timeSlider);
		this.domElement.append(processHtml);
		timeSlider.bind('change', {
			processGroup: processGroup
		}, this.handleTimeChange);

		if (processGroup.length > 1) {
			var processSlider = $(processSliderTemplate);
			timeSlider.after(processSlider);
			processSlider.bind('change', {
				processGroup: processGroup,
				scope:this
			}, this.handleProcessChange);

			for (var i = 0; i< processGroup.length; ++i) {
				// Modify slider / drop down for each process
			}
			// bind event to second slider
		}
	}
}
ProcessToolbox.prototype.handleTimeChange = function(event) {
	var dataId = $(this).val();
	var processGroup = event.data.processGroup;
	var process = processGroup[0];
	process.render(dataId);
}
ProcessToolbox.prototype.handleProcessChange = function(event) {
	var processId = $(this).val();
	var processGroup = event.data.processGroup;
	scope.updateSliderFromProcess(processGroup[processId], $(this));
}
ProcessToolbox.prototype.updateSliderFromProcess = function(process, slider) {
	var alpha;
}