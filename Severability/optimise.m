%inst=3;

max_size = 60;

for inst=1:5
for graph_size=[5000]
for mixing=[0.1:0.1:0.6]
for time=[2:4]
	%if ((inst==1)&&(time==2))
	%	break
	%end

	disp({'inst:',inst,'size:',graph_size,'mixing:', mixing,'time', time})
	% Make sure the max_size value makes sense.
	if ((max_size>graph_size) || (max_size<0))
		max_size = graph_size;
	end
	% Other string filenames
	network_file = strcat('benchmark_2_1/',num2str(graph_size),'_',num2str(inst),'/',num2str(mixing),'/network.dat');
	realcomm_file = strcat('benchmark_2_1/',num2str(graph_size),'_',num2str(inst),'/',num2str(mixing),'/community.dat');
	outputname = strcat('sev',num2str(graph_size),'_',num2str(inst),'m',num2str(mixing),'t',num2str(time),'r',num2str(max_size),'.mat');

	adj_matrix=adjacency_matrix(network_file);
	realcomm=load(realcomm_file);

	[ commrun node commsize sev output commrun_small separated missing overlap ...
	    foundcomm] = sev_multi_level( adj_matrix, realcomm, time, max_size );

	disp({'singletons:',length(missing),'overlap',length(overlap)})
	mutualinformation = norm_mutualinfo(realcomm(:,2),foundcomm)
	outputline={'inst:',inst,'size:',graph_size,'mixing:', mixing,'time', time,'singletons:',length(missing),'overlap',length(overlap),'mutualinformation',mutualinformation}
	save(outputname)
	if (length(missing)+length(overlap))==0
		break
	end
	%save sev500m03t1.mat
end
end
end
end
