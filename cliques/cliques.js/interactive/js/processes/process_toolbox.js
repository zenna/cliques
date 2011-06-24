dynamicProcessTemplate = "<div><span class='processName'></span>\
<input class = 'timeSlider' name='r' type='range' min='0' max='529' value='0'>\
<span class='time'>a</span>\
</div>\
<br/>\
"

processSliderTemplate = "<input class = 'processSlider' name='r' type='range' min='0' max='10' value='0'/>"

var ProcessToolbox = function(view) {
	this.domElement = $("<div></div>").addClass('toolbox');
	$(view.renderer.domElement).parent().append(this.domElement);
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
	var process = processGroup[0];
	console.log(dataId);
	var time = process.data[dataId]['time'];
	$(this).next('.time').text(time);
	
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

louvain = {
	type:'louvain',
	dependencies: {
		type:'stability',
		time:0.234
	},
	data:[{
		partition:'034202',
		action:'consider',
		time:0
	},{
		partition:'034202',
		action:'consider',
		time:0
	},{
		partition:'034202',
		action:'position',
		time:1
	},
	]
}

louvain2 = {
	type:'louvain',
	dependencies: {
		type:'stability',
		time:4.3
	},
	data:[{
		partition:'034202',
		action:'consider',
		time:0
	},{
		partition:'034202',
		action:'consider',
		time:0
	},{
		partition:'034202',
		action:'position',
		time:1
	},
	]
}
