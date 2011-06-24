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
	this.graph.move_nodes(norm_values, 0, 1);
	this.graph.match_edge_colours_to_node();


	// var alpha = {
	// type:'basin',
	// dependencies:['stability'],
	// data:[{time:0.234, values:[{basin:0,nodes:[0,2,3]},
	// {basin:1,nodes:[3,4,5]}];
	// this.graph.match_edge_colours_to_node();
	// }

}