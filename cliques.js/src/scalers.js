var cliques = cliques || {};

cliques.Scaler = function(minValue, maxValue) {
	this.minValue = typeof(minValue) != 'undefined'  ? minValue : Number.MAX_VALUE;
	this.maxValue = typeof(minValue) != 'undefined'  ? minValue : Number.MIN_VALUE;
}
cliques.Scaler.prototype.updateExtrema = function(updateValue) {
	this.updateMaxima(updateValue);
	this.updateMinima(updateValue);
}
cliques.Scaler.prototype.updateMaxima = function(updateValue) {
	if (updateValue > this.maxValue) {
		this.maxValue = updateValue;
		doUpdateRanges = true;
	}
}
cliques.Scaler.prototype.updateMinima = function(updateValue) {
	if (updateValue < this.minValue) {
		this.minValue = updateValue;
		doUpdateRanges = true;
	}
}
cliques.Scaler.prototype.updateRanges = function() {
	this.shift = 0.0 - this.minValue;
	this.range = this.maxValue - this.minValue;
	this.doUpdateRanges = false;
}
cliques.LinearScaler = function() {
	cliques.Scaler.call(this);
}
cliques.LinearScaler.prototype = new cliques.Scaler();
cliques.LinearScaler.prototype.constructor = cliques.LinearScaler;

cliques.LinearScaler.prototype.scaleValue = function(value) {
	if (this.doUpdateRanges = true) {
		this.updateRanges();
	}
	return (value + this.shift)/this.range;
}