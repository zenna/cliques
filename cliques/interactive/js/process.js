stabilities = {
	type:'stability',
	data:[{
		time:0.1,
		values:[1.3,-2.3,-0.2]
	},{
		time:0.6,
		values:[9.3,-0.3,-1.2]
	}]
};

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

var NodeProcess = function(process, graph) {
	this.type = process.type;
	this.data = process.data;
	this.graph = graph;
};
NodeProcess.prototype.render = function(time, data_id) {
	//WICKED YOU KNOW!!!
	alert("rendering"+data_id);
}
NodeProcess.prototype.hide = function(time) {

}
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
staticProcessTemplate = "<div><span class='processName'></span>@\
<input class = 'time' name='r' type='range' min='0' max='10' value='0'>\
<span class='value'></span>\
</div>\
"

dynamicProcessTemplate = "<div><span class='processName'></span>@\
<select class='processes'></select>\
<input class = 'time' name='r' type='range' min='0' max='10' value='0'>\
<span class='value'></span>\
</div>\
"

ProcessToolbox.prototype.updateDomElement = function() {
	// Generate html from processes

	for (var process in this.processes) {
		var processHtml = $(dynamicProcessTemplate);
		// for each process group
		//create a template
		// for

		this.domElement.html('');
		for (var i=0;i<this.processes.length;++i) {
			var process = this.processes[i];
			var processHtml = $(dynamicProcessTemplate);
			processHtml.children('.processName').text(process.type);

			if (typeof process.dependecies != 'undefined') {
				for (var j=0;j<process.data.length;++j) {
					var current_data = process.data[j];
					var option = $("<option value=''></option>");
					option.text(current_data.time);
					option.val(j);
					processHtml.children('.processes').append(option);
				}
			}
			this.domElement.append(processHtml);
			processHtml.bind('change', {
				process:process
			}, this.handleChange);
		}
	}
}
ProcessToolbox.prototype.handleChange = function(event) {
	var data_id = $(this).find('select option:selected').val();
	var time = $(this).find('.time').val();
	var process = event.data.process;
	process.render(time,data_id);
}