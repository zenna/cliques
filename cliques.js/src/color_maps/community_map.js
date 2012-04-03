cliques.CommunityColorMap = function() {
	cliques.ColorMap.call(this);
}
cliques.CommunityColorMap.prototype = new cliques.ColorMap;
cliques.CommunityColorMap.prototype.constructor = cliques.CommunityColorMap;

cliques.CommunityColorMap.prototype.getColors = function(community) {
	var groups = [];
	var colors = [];
	this.scaler.updateExtrema(0);
	this.scaler.updateExtrema(1);
	
	for (var i=0;i<community.length;++i) {
		var scaledColor = this.scaler.scaleValue(community[i]);
		var color = this.mapColor(scaledColor);
		colors[i] = color;
	}
	return colors;
}