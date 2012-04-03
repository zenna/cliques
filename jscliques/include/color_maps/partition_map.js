cliques.PartitionColorMap = function() {
	cliques.ColorMap.call(this);
}
cliques.PartitionColorMap.prototype = new cliques.ColorMap;
cliques.PartitionColorMap.prototype.constructor = cliques.PartitionColorMap;

cliques.PartitionColorMap.prototype.getColors = function(partition) {
	var groups = [];
	var colors = [];
	this.scaler.updateExtrema(0);
	//this.scaler.updateExtrema(partition.length);
	for (var i =0;i<partition.length;++i) {
		if (!(groups[partition[i]] instanceof Array)) {
			groups[partition[i]] = [];
		}
		groups[partition[i]].push(i);
	}
	this.scaler.updateExtrema(groups.length -1);
	groups.sort(sortBySize)

	function sortBySize(a,b) {
		return b.length - a.length;
	}
	
	for (var i =0;i<groups.length; ++i) {
		var scaled_color = this.scaler.scaleValue(i);
		var color = this.mapColor(scaled_color);
		for (var j=0;j<groups[i].length; ++j) {
			var node = groups[i][j];
			colors[node] = color;
		}
	}
	return colors;
}