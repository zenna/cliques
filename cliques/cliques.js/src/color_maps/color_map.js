// Need colormap to:
// a) given partition, return list of colours corresponding to each node
// b) given a partition, modify materials object to colour items

cliques.ColorMap = function() {
	this.scaler = new cliques.LinearScaler();
	this.mapType = 'michael';
	// string argument vs inherited
	// string means will need case statement but can switch type easily
}

cliques.ColorMap.prototype.updateExtrema = function(value) {
	this.scaler.updateExtrema(value);
}
cliques.ColorMap.prototype.render = function() {
	// Render this colourmap to html object
}

cliques.ColorMap.prototype.mapColor = function(value) {
	switch(this.mapType) {
		case 'grey':
		return [value, value, value]
		break;
		case 'michael':
		var r,g,b;
		if (value < 0.25) {
			b = 1.0;
			g = 4*value;
			r = 0.0;
		} else if (value < 0.5) {
			b = 2 - 4 * value;
			g = 1.0;
			r = 0.0;
		} else if (value < 0.75) {
			g = 1.0;
			r = 4 * value - 2;
			b = 0.0;
		}
		else {
			r = 1.0;
			g = 4 - 4*value;
			b = 0.0;
		}
		return [r, g, b];
		break;
		case 'reds':
		var r,g,b;
		if (value < 0.5) {
			b = 1- 2*value;
			g = 2*value;
			r = 0.0;
		} else {
			g = 2 - 2 * value;
			r = 2* value -1 ;
			b = 0.0
		}
		return [r, g, b];
		break;
		case 'hsv':
		return cliques.hsvToRgb([value, 1.0, 1.0]);
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

cliques.ColorMap.prototype.scaleColor = function(color, scaleFactor) {
	var scaledColor = [];
	scaledColor[0] = color[0] * scaleFactor;
	scaledColor[1] = color[1] * scaleFactor;
	scaledColor[2] = color[2] * scaleFactor;
	return scaledColor;
}

cliques.ColorMap.prototype.add = function(color1, color2) {	
	var combinedColor = [];
	combinedColor[0] = color1[0] + color2[0];
	combinedColor[1] = color1[1] + color2[1];
	combinedColor[2] = color1[2] + color2[2];
	return combinedColor;
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

cliques.hsvToRgb = function(hsv) {
	var h = hsv[0] * 360;
	var s = hsv[1];
	var v = hsv[2]
	// var chroma = hsv[2] * hsv[1];
	// var h_d = hsv[0] * 6.0;
	  // import math
    var hi = Math.floor(h / 60.0) % 6;
    var f =  (h / 60.0) - Math.floor(h / 60.0)
    var p = v * (1.0 - s)
    var q = v * (1.0 - (f*s))
    var t = v * (1.0 - ((1.0 - f) * s))
    var rgb = {
        0: [v, t, p],
        1: [q, v, p],
        2: [p, v, t],
        3: [p, q, v],
        4: [t, p, v],
        5: [v, p, q]
    };
    return rgb[hi];	
}
