cliques.EnergyColorMap = function() {
	cliques.ColorMap.call(this);
}
cliques.EnergyColorMap.prototype = new cliques.ColorMap;
cliques.EnergyColorMap.prototype.constructor = cliques.EnergyColorMap;

cliques.EnergyColorMap.prototype.getColor = function(value, scaleFactor) {
	var color = this.mapColor(this.scaler.scaleValue(value));
	if (typeof(scaleFactor) == 'undefined') {
		return color;
	}
	else {
		return this.scaleColor(color);
	}
}