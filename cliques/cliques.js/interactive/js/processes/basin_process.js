basinsMetaTemplate = "<div class ='meta'>time:<span class='time'></span>\
Num Basins:<span class = 'numBasins'></span>\
<ul class = 'basinOptions'></ul>\
</div>"

basinOptionsTemplate = "<li>Maxima Id:<span class='maximaId'></span>\
Basin Size:<span class ='basinSize'></span><input type='checkbox'>isolate</input>, focus</li>"

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
		for (var j = 0;j<values[i].nodes.length;++j) {
			var node_id = values[i].nodes[j];
			var probability = values[i].metas[j];
			var scaledBasinColor = colorMap.scaleColor(basinColor,probability);
			if (node_id == 5133) {
				console.log("found you");
			}
			if (typeof(rgbs[node_id]) == 'undefined') {
				rgbs[node_id] = [0.0,0.0,0.0];
			}
			var combinedColor = colorMap.add(rgbs[node_id], scaledBasinColor);
			rgbs[node_id] = combinedColor;
			norm_values[node_id] = colorMap.scaler.scaleValue(i);
		}
		// graph.highlight_node(values['basin'], basinColor));
		rgbs[values[i]['basin']] = [1.0,1.0,1.0];

	}
	
	this.graph.paint_nodes(rgbs);
	//this.graph.move_nodes(norm_values, 0, 1);
	this.graph.match_edge_colours_to_node();

	// var alpha = {
	// type:'basin',
	// dependencies:['stability'],
	// data:[{time:0.234, values:[{basin:0,nodes:[0,2,3]},
	// {basin:1,nodes:[3,4,5]}];
	// this.graph.match_edge_colours_to_node();
	// }

}

BasinProcess.prototype.makeHtmlFromBasin = function(maximaId, basinSize) {
	var basinHtml = $(basinOptionsTemplate);
	basinHtml.attr('id',maximaId);
	basinHtml.find('.maximaId').text(maximaId);
	basinHtml.find('.basinSize').text(basinSize);
	

	return basinHtml;
}

BasinProcess.prototype.updateMeta = function(dataId) {
	var data = this.data[dataId];
	var time = data['time'];
	var numBasins = data.values.length;
	
	var basinsMetaHtml = $(basinsMetaTemplate);
	basinsMetaHtml.find('.time').html(time);
	basinsMetaHtml.find('.numBasins').html(numBasins);
	var basinOptions = basinsMetaHtml.find('.basinOptions');
	
	for (var i=0;i<data.values.length;++i) {
		var maximaId = data.values[i]['basin'];
		var basinSize = data.values[i]['nodes'].length;
		var basinHtml = this.makeHtmlFromBasin(maximaId, basinSize);
		basinOptions.append(basinHtml);
	}
	var self = this;
	basinOptions.bind('change', {
		scope: this, data:data, dataId:dataId
	}, function(event) {self.handleIsolateBasin(event)});
		
	return basinsMetaHtml;
}

BasinProcess.prototype.handleIsolateBasin = function(event) {
	var self = this;
	var data = event.data.data;
	this.render(event.data.dataId);
	var rgbs = [];
	var colors = event.data.scope.graph.nodes.colors;
	for (var i=0;i<colors.length;++i) {
		rgbs.push([0.0,0.0,0.0]);
	}
	var checkboxes = event.currentTarget;
	var checkedBoxes = $(checkboxes).find(':checked');
	checkedBoxes.each(function(index, element) {
		var maximaId = parseInt($(element).parent().attr('id'));
		var maximaIndex;
		for (var i=0;i<data.values.length;++i) {
			if (data.values[i]['basin'] == maximaId) {
				maximaIndex = i;
				break;
			}
		}
		for (var i=0;i<data.values[maximaIndex]['nodes'].length;++i) {
			var node = data.values[maximaIndex]['nodes'][i];
			rgbs[node] = [colors[node].r, colors[node].g, colors[node].b];
		}
	})
	
	
	if (checkedBoxes.length > 0) {
		this.graph.paint_nodes(rgbs);
		this.graph.match_edge_colours_to_node();
	}
}

BasinProcess.prototype.showNodeData = function(nodeId) {
	return "I'm a basin process stabilty of node" + nodeId;
}
