dynamicProcessTemplate = "<div><span class='processName'></span>\
<input class = 'timeSlider' name='r' type='range' min='0' max='529' value='0'>\
<span class='time'></span>\
<div class = 'meta'></div>\
</div>\
<br/>\
";

processSliderTemplate = "<input class = 'processSlider' name='r' type='range' min='0' max='10' value='0'/>";

processPopup = "<div class='processPopup'></div>";

var ProcessToolbox = function(view) {
	this.domElement = $("<div></div>").addClass('toolbox');
	$(view.renderer.domElement).parent().append(this.domElement);
	this.processPopup = $("<div></div>").addClass('processPopup');
	$(view.renderer.domElement).parent().append(this.processPopup);
	this.processes = {};
}
ProcessToolbox.prototype.addProcess = function(process) {
	var type = process.type;
	if (typeof this.processes[type] == "undefined") {
		this.processes[type] = [];
	}
	this.processes[type].push(process);
	this.updateDomElement();
}

ProcessToolbox.prototype.updateDomElement = function() {
	// Generate html from processes

	this.domElement.html('');
	for (var process in this.processes) {
		var processGroup = this.processes[process];
		var processHtml = $(dynamicProcessTemplate);
		processHtml.children('.processName').text(process);
		var timeSlider = processHtml.children('input.timeSlider')
		this.updateSliderFromProcess(processGroup[0], timeSlider);
		this.domElement.append(processHtml);
		timeSlider.bind('change', {
			processGroup: processGroup
		}, this.handleTimeChange);

		if (processGroup.length > 1) {
			var processSlider = $(processSliderTemplate);
			timeSlider.after(processSlider);
			processSlider.bind('change', {
				processGroup: processGroup,
				scope:this
			}, this.handleProcessChange);

			for (var i = 0; i< processGroup.length; ++i) {
				// Modify slider / drop down for each process
			}
			// bind event to second slider
		}
	}
}
ProcessToolbox.prototype.handleTimeChange = function(event) {
	var dataId = $(this).val();
	var processGroup = event.data.processGroup;
	var process = processGroup[0]; // This will fail for a dynamic process dependent on stability
	var time = process.data[dataId]['time'];
	process['currentTime'] = dataId;
	$(this).next('.time').text(time);
	var metaHtml = process.updateMeta(dataId);
	$(this).siblings('.meta').html(metaHtml);
	process.render(dataId);
}
ProcessToolbox.prototype.handleProcessChange = function(event) {
	var processId = $(this).val();
	var processGroup = event.data.processGroup;
	scope.updateSliderFromProcess(processGroup[processId], $(this));
}
ProcessToolbox.prototype.updateSliderFromProcess = function(process, slider) {
	var num_steps = process.data.length - 1;
	slider.attr('max', num_steps);
}

ProcessToolbox.prototype.popupNodeDataFunctions = function(data) {
	// Call
	var nodeDataHtml = "";
	
	for (var i=0; i<this.processes['stability'].length;++i) {
		var currentTime = this.processes['stability'][i]['currentTime'];
		nodeDataHtml += this.processes['stability'][i].showNodeData(data,currentTime);
	}
	this.processPopup.text(nodeDataHtml);
}