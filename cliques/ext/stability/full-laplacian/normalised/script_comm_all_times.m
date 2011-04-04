% Inputs: Graph saved in data.mat in the form: 
% 		"source1 dest1 weight1" 
% 		"source2 dest2 weight2" 
% 		...
%
% Ouptut: files starting with the string defined in name_out
%		+ *.dat => partitioning (node partition)
%		+ *.mat => distribution of the variation of informations
%		+ *.stdout => text file containing in columns:
%			1) Markov time	
%			2) Best value of stability
%			3) Mean stability
%			4) Std deviation of stability
%			5) # communities for best partition
%			6) Mean # communities	
%			7) Std deviation # communities
%			8) Mean VI 	
%			9) Std deviation of VI 	
%			10) Median VI	
%			11) Min VI		
%			12) Max VI
%

timesx = []
% Preprocessing3um
% 
% folder_root = '/home/zenna/repos/uniproj/data/brains30/'
% basename = 'COSNV'
% out_folder = strcat('/home/zenna/repos/uniproj/data/brain_new_stabilities/',basename,num2str(iii))
% filename = strcat(folder_root,basename,num2str(iii),'.txt')

% timesx = []
% % Preprocessing3um
%
%iii=4
folder_root = '/home/zenna/repos/uniproj/data/brains30/'
basename = 'COSNV'
out_folder = strcat('/home/zenna/repos/uniproj/data/random_stabilities/','R',basename,num2str(iii))
%out_folder = '/home/zenna/repos/uniproj/data/graph_stabilities/random_brain2'
filename = strcat(folder_root,basename,num2str(iii),'.txt')
%filename = '/home/zenna/repos/uniproj/data/graphs/random_n12.mlab'

Graph = load(filename)
%Graph = COSNVRANDOM%load(COSNVRANDOM);

%Randomise brain
for i = 1:10000000
    n = floor(rand * 9900) + 1;
    p = floor(rand * 9900) + 1;
    buffer = Graph(n,3);
    Graph(n,3) = Graph(p,3);
    Graph(p,3) = buffer;
end

%Graph = Graph2;
%Graph = new_brain;
nb_nodes = max(Graph(:,1))+1;
l=size(Graph,1);
Graph(:,1) = Graph(:,1)+1;
Graph(:,2) = Graph(:,2)+1;
m=max(max(Graph(:,1)),max(Graph(:,2)));
A=sparse(m,m);
for i=1:l
   A(Graph(i,1),Graph(i,2))=Graph(i,3);
end
clear Graph;
A=A+A';
save A.mat A;
clear A;


% Parameters
name_out = out_folder; % Name of output files
time = 0.008; % Starting time
nTrials = 100; % Number of restarts for louvain algorithm
distr_freq = Inf; % Frequency at which the distribution of the variation of informations is computed. ex: 10 => every 10 samples

% Initialisation
stability_list=zeros(nTrials,1);
nb_comm_list=zeros(nTrials,1);
condition=true;
compteur=0;
mkdir(name_out)


%Loop on Markov times
while(condition)
    compteur=compteur+1;
    lnk = zeros(nb_nodes, nTrials);


    disp('Exponential computation');
    tstart = tic;
    graph = expmat_revised(time);
    toc(tstart)


    disp('Partitioning');
    tstart = tic;
    stability_best = -Inf;
    communities_best = [];
    nb_comm_best = Inf;
        
    for i=1:nTrials
        [stability, nb_comm, communities] = main_community_mex(graph);
        if nb_comm == Inf
            nb_comm=1;
            condition=false;
            break;
        end
        nb_comm_list(i)=nb_comm;
        stability_list(i)=stability;
        communities = main_hierarchy(communities);
        lnk(:,i) = communities(:,2);
        if stability>stability_best
            stability_best = stability;
            communities_best = communities;
            nb_comm_best = nb_comm;
        end
    end
    toc(tstart)
    
    communities = communities_best;
    clear communities_best;

    disp('Variation of Information');
    if mod(compteur,distr_freq)==0
    	[var_info, var_vi, median_vi, min_vi, max_vi, matrix, distr] = main_varinfo_new(lnk);
	clear matrix;
	save(['distr_' num2str(time) '.mat'],'distr');
    else
	[var_info, var_vi, median_vi, min_vi, max_vi] = main_varinfo_new(lnk);
    end
    clear lnk;
    
    stability_best
    nb_comm_best
    %var_info
    
    name=[name_out '/' num2str(time,'%10.8f')];
    out_id=fopen(name,'w');
    for i=1:length(communities(:,1))
        fprintf(out_id,int2str(communities(i,1)));
        fprintf(out_id,'\t');
        fprintf(out_id,int2str(communities(i,2)));
        fprintf(out_id,'\n');
    end
    fclose(out_id);
    
    name=[name_out '.stdout'];
    out_id=fopen(name,'a');
    fprintf(out_id,num2str(time,'%10.8f'));
    fprintf(out_id,'\t');
    fprintf(out_id,num2str(stability_best,'%10.8f'));
    fprintf(out_id,'\t');
    fprintf(out_id,num2str(mean(stability_list),'%10.8f'));
    fprintf(out_id,'\t');
    fprintf(out_id,num2str(std(stability_list),'%10.8f'));
    fprintf(out_id,'\t');
    fprintf(out_id,int2str(nb_comm_best));
    fprintf(out_id,'\t');
    fprintf(out_id,num2str(mean(nb_comm_list),'%10.8f'));    
    fprintf(out_id,'\t');
    fprintf(out_id,num2str(std(nb_comm_list),'%10.8f'));  
    fprintf(out_id,'\t');
    fprintf(out_id,num2str(var_info));
    fprintf(out_id,'\t');
    fprintf(out_id,num2str(var_vi));
    fprintf(out_id,'\t');
    fprintf(out_id,num2str(median_vi));
    fprintf(out_id,'\t');
    fprintf(out_id,num2str(min_vi));
    fprintf(out_id,'\t');
    fprintf(out_id,num2str(max_vi));
    fprintf(out_id,'\n');
    fclose(out_id);
    
    %timesx = [timesx time];

    step = time/50;
    time = time+time/25
    
    if(time>10000)
        condition=false;
    end
end
fclose('all');
if doExit == 1
 exit
end
