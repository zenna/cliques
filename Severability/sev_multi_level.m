function [ commrun node commsize sev output commrun_small separated missing overlap ...
    foundcomm] = sev_multi_level( adj_matrix, realcomm, time, radius )
%SEV_MULTI_LEVEL Iterates sev_one_subcomm to optimise severability
%   Very similar in spirit to the Louvain aggregation algorithm for
%   modularity. Calls "sev_one_subcomm" repeatedly to get a covering of the
%   graph by communities. Then collapses all of the communities, and
%   repeats.
%
%   adj_matrix: symmetric n*n matrix detailing connectivity of graph
%   time:   can be singleton
%   realcomm:   community.dat file specifying the "actual" communities

%% Main "simulation"
graph_size=length(adj_matrix);
commrun=zeros(graph_size,graph_size); % Contains all the communities found

covered=zeros(1,graph_size);  % List of all the nodes that have already been covered
node=zeros(1,graph_size);    % Node associated with "i"th entry
commsize=zeros(1,graph_size);   % Sizes of the communities (redundant with commrun)
sev=zeros(1,graph_size);    % Severabilities of the communities found

curr_entry = 0;
for i=1:graph_size;
    if covered(i)==0
        curr_entry=curr_entry+1;
        node(curr_entry)=i;
        [comm commsize(curr_entry) sev(curr_entry)]=sev_node(i,adj_matrix,time,radius);

    % Put entries in covered
    if commsize(curr_entry)~=graph_size
        for j=1:(commsize(curr_entry))
            covered(comm(j))=covered(comm(j))+1;
        end
    end
    
    %comm
    commrun(curr_entry,:)=comm;
    fprintf('\t: %4d %4d %4d %4f \n', ...
        [sum((covered~=0)) i commsize(curr_entry) sev(curr_entry)]);

    end    
end;
% Resize all of the arrays
commrun=commrun(1:curr_entry,:);
node=node(1:curr_entry);
commsize=commsize(1:curr_entry);
sev=sev(1:curr_entry);

%% Data Analysis
output=zeros(length(node),5);
for i=1:length(node)
    comm=sort(commrun(i,1:commsize(i)));
    
    % Gets the list of communities that each of the elemnts in "comm"
    % belong to. (There should be a lot of repeats).
    if var(realcomm(comm,2))~=0
        if commsize(i)==graph_size
            output(i,:)=-2; % Entire graph
        else
            output(i,:)=-1; % Mixed communities
        end
    else
        real_community=find(realcomm(:,2)==mean(realcomm(comm,2)));
        difference=setdiff(real_community,comm);
        output(i,1:length(difference))=difference;
    end
end
separated=unique(output(output>0));

% Get rid of all communities that are the whole graph
commrun_small=commrun;
commrun_small(commsize==graph_size,:)=[];

%% Partition generation
% Create a partition from the communities found
%graph_size=1000;
foundcomm = zeros(1,graph_size);
for i=1:graph_size
    if (mod(i,10)==0)
        if i~=10
            fprintf('\b\b\b\b')
        else
            fprintf('.')
        end
        fprintf('%4d',i)    % Progress bar
    end
    for j=1:size(commrun_small,1)
        if ismember(i,commrun_small(j,:))
            foundcomm(i)=j;
            break
        end
        foundcomm(i)=-i;
    end
end
%%

% Find how well the communities cover the space.
% Also find which nodes are repeated in multiple communities:

% First make commrun_small smaller again, taking out all communities that aren't used at all
commrun_small=commrun_small(unique(foundcomm(foundcomm>0)),:);

missing=[];
overlap=[];
for i=1:graph_size
    if isempty(find(commrun_small==i,1,'first'))
        missing(end+1)=i;
    end
    if length(commrun_small(commrun_small==i))>1
        overlap(end+1)=i;
    end
end

%%
end

