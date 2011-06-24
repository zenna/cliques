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