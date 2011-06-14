{
	process: [{
		type:'louvain',
		dependencies:{type:'stability', time:0.234},
		data:[{partition:'034202',action:'consider', time:0},
			{partition:'034202',action:'consider', time:0},
			{partition:'034202',action:'position', time:1},
		]
	},
	{
		type:'k_lin',
		dependencies:{type:'stability', time:0.234},
		data:[{partition:'034202',action:'consider', time:0},
			{partition:'034202',action:'consider', time:0},
			{partition:'034202',action:'position', time:0},
		]
	},
	{
		type:'stability',
		data:[{time:0.234,values:[1.3,-2.3,-0.2]},
		{time:0.239,values:[9.3,-0.3,-1.2]}]
	},
	{
		type:'maxima',
		dependencies:{type:'stability', time:0.234},
		data:[{nodes:[0,1,2,3,4,10]}]
	},
	{
		type:'bottlenecks',
		dependencies:{type:'stability', time:0.234},
		data:[{nodes:[0,1,2,3,4,10]}]
	},
	{
		type:'basin',
		dependencies:{type:'stability', time:0.234},
		data:[{basin:0,nodes:[0,2,3]},
		{basin:1,nodes:[3,4,5]}]
	},
	{
		type:'pbasin',
		dependencies:{type:'stability', time:0.234},
		data:[{basin:0,nodes:[0,2,3]},
		{basin:1,nodes:[3,4,5]}]
	},
	]
}

// Generate data
// Write the parsers
// Create the UI



//id vs partititon
// id: smaller data size, simple lookup
// partition: can move process another landscape, e.g. sampled (what use is this?)
// not susceptible to order in which partitions were generated. Possible sol: ensure processes are only made from partition data
