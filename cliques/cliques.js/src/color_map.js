cliques.Scaler = function(minValue, maxValue) {
    this.minValue = typeof(minValue) != 'undefined'  ? minValue : Number.MAX_VALUE;
    this.maxValue = typeof(minValue) != 'undefined'  ? minValue : Number.MIN_VALUE;
}

cliques.Scaler.prototype.updateExtreme = function(updateValue) {
    if (update_value > this.maxValue) {
        this.maxValue = update_value;
        doUpdateRanges = true;
    }
}

cliques.Scaler.prototype.updateRanges = function() {
    this.shift = 0.0 - this.minValue;
    this.range = this.maxValue - this.minValue;
    this.doUpdateRanges = false;
}

cliques.LinearScaler = function() {
    this.prototype = cliques.Scaler;
}

cliques.LinearScaler.prototype.scale_value() = function(value) {
    if (this.doUpdateRanges = true) {
        this.updateRanges();
    }
    return (value + this.shift)/this.range;
}

cliques.ColorMap = function() {
    this.scaler = new scaler()
    // string argument vs inherited
    // string means will need case statement but can switch type easily
}

cliques.ColorMap.prototype.render = function() {
	// Render this colourmap to html object
}

cliques.ColorMap.prototype.getColorRGB = function() {
	// find scaled value
	// convert to colour as RGB
}

cliques.ColorMap.prototype.getColorHex = function() {
	//get colour as RGB
	//cliques.RGBToHex
}

cliques.convertColorSpace = function() {
	
}

cliques.hexToRGB = function() {
	
}

cliques.HSVToRGB = function() {
	
}
