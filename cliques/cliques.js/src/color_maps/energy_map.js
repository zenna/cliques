cliques.EnergyColorMap = function() {
	cliques.ColorMap.call(this);
}
cliques.EnergyColorMap.prototype = new cliques.ColorMap;
cliques.EnergyColorMap.prototype.constructor = cliques.EnergyColorMap;

cliques.EnergyColorMap.prototype.getColor = function(value) {
	return this.mapColor(this.scaler.scaleValue(value));
}