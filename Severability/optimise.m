
push_ahead=-1;
radius = 50;

for graph_size=[500]
for mixing=0.5
for inst=5
for time=[2]

%	if ((mixing==0.1&&inst<5))
%		break
%	end

	disp({'inst:',inst,'size:',graph_size,'mixing:', mixing,'time', time})
	% Make sure the radius value makes sense.
	if ((radius>graph_size) || (radius<0))
		radius = graph_size;
	end
	% Other string filenames
	network_file = strcat('benchmark_2_1/',num2str(graph_size),'_',num2str(inst),'/',num2str(mixing),'/network.dat');
	realcomm_file = strcat('benchmark_2_1/',num2str(graph_size),'_',num2str(inst),'/',num2str(mixing),'/community.dat');
	outputname = strcat('sev',num2str(graph_size),'_',num2str(inst),'m',num2str(mixing),'t',num2str(time),'r',num2str(radius),'.mat');

	adj_matrix=adjacency_matrix(network_file);
	realcomm=load(realcomm_file);

	[ commrun node commsize sev output commrun_small separated missing overlap ...
	    foundcomm] = sev_multi_level( adj_matrix, realcomm, time, push_ahead, radius );

	disp({'singletons:',length(missing),'overlap',length(overlap)})
	mutualinformation = norm_mutualinfo(realcomm(:,2),foundcomm)
	save(outputname)
	if (length(missing)+length(overlap))==0
		break
	end
	%save sev500m03t1.mat
end
end
end
end
