function [ commrun commsize sev commrun_small ...
    missing overlap foundcomm] = ...
    sev_partition( adj_matrix, realcomm, time, max_size )
%SEV_PARTITION Softly "partitions" a graph using severability
%   [commrun commsize sev commrun_small ...
%      missing overlap foundcomm] = ...
%   sev_partition( adj_matrix, realcomm, time, max_size)
%
%   adj_matrix: adjacency matrix
%     realcomm: community.dat file specifying the "actual" communities
%         time: Markov time
%     max_size: size limit in the community of the search algorithm
%

%      commrun: k by [graph_size] matrix, specifying the communities
%               commrun(i,1:commsize(i)) gives you all community members
%               in community i.
%     commsize: the sizes of the community found (length k vector)
%          sev: severability of the found communities (length k vector)
%commrun_small: Reduced commrun, taking out communities unused in the
%               partitioning.
%      missing: All the nodes that do not belong to a found community (orphans)
%      overlap: All the nodes that belong to more than one community in the list
%    foundcomm: The partition scheme.
%   
%   Copyright (c) 2010-2011 Yun William Yu
%   Revision 2011-06-22


%% Main "simulation"
graph_size=length(adj_matrix);
commrun=zeros(graph_size,graph_size); % Contains all the communities found

covered=zeros(1,graph_size);  % List of all the nodes that have already been covered
node=zeros(1,graph_size);    % Node associated with "i"th entry
commsize=zeros(1,graph_size);   % Sizes of the communities (redundant with commrun)
sev=zeros(1,graph_size);    % Severabilities of the communities found

curr_entry = 0;
for i=1:graph_size;

    % if the node i is on average more connected to the already covered, then skip.
    edges_out=adj_matrix(i,:);
    if sum(edges_out(covered==0))<sum(edges_out(covered~=0))
        forceorphan=true;
    else
        forceorphan=false;
    end
    if (covered(i)==0) && (forceorphan==false)
        curr_entry=curr_entry+1;
        node(curr_entry)=i;
        [comm commsize(curr_entry) sev(curr_entry)]=sev_node(i,adj_matrix,time,max_size);

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

