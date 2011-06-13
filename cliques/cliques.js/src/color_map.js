// Need colormap to:
// a) given partition, return list of colours corresponding to each node
// b) given a partition, modify materials object to colour items

cliques.ColorMap = function() {
	this.scaler = new cliques.LinearScaler();
	this.mapType = 'grey';
	// string argument vs inherited
	// string means will need case statement but can switch type easily
}
cliques.ColorMap.prototype.render = function() {
	// Render this colourmap to html object
}

cliques.ColorMap.prototype.mapColor = function(value) {
	switch(this.mapType) {
		case 'grey':
		return [value, value, value]
		break;
		case 'reds':
		return [value, 1.0, 1.0];
	}
}

cliques.ColorMap.prototype.getColorRGB = function() {
	// find scaled value
	// convert to colour as RGB
}
cliques.ColorMap.prototype.getColorHex = function() {
	//get colour as RGB
	//cliques.RGBToHex
}
cliques.PartititionColorMap = function() {
	cliques.ColorMap.call(this);
}
cliques.PartititionColorMap.prototype = new cliques.ColorMap;
cliques.PartititionColorMap.prototype.constructor = cliques.PartititionColorMap;

cliques.PartititionColorMap.prototype.getColors = function(partition) {
	var groups = [];
	var colors = [];
	this.scaler.updateExtrema(0);
	this.scaler.updateExtrema(partition.length);
	for (var i =0;i<partition.length;++i) {
		if (!(groups[partition[i]] instanceof Array)) {
			groups[partition[i]] = [];
		}
		groups[partition[i]].push(i);
	}
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
cliques.convertColorSpace = function() {

}
cliques.hexToRGB = function() {

}
cliques.HSVToRGB = function() {

}
cliques.intToRgb = function(integer) {
	var n = Math.floor(integer / 256);
	var r = integer % 256;
	var n2 = Math.floor(n / 256);
	var r2 = n % 256;
	return [n2,r2,r];
}
cliques.rgbToInt = function(rgb) {
	return rgb[0] * 256 * 256 + rgb[1] * 256 + rgb[2];
}