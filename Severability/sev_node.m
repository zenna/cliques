function [ commrun commsize sev ] = sev_node( start_node, adj_matrix, t, max_size)
%SEV_NODE Greedy optimisation for the sev-community of a node, with EM
%   Will always give you a community with the original node in it.
%   Copyright (c) 2010 Yun William Yu
commrun=[];
% Doesn't make sense for the max size to be greater than the community size
max_size = min(max_size,size(adj_matrix,1));

%%
graph_size=length(adj_matrix);
%A = adj_matrix(community,community);
D = diag(sum(adj_matrix)); % diagonal degree matrix
P = (D^-1)*adj_matrix;  % Markovian transition matrix
%P = full(P);

% Get neighbors of startnode
adj_nodes=find(sum(adj_matrix(start_node,:),1));
neighbors=setdiff(adj_nodes,start_node);
sev_orders = zeros(length(neighbors),1);      
for i=1:length(neighbors)
    Q = P([start_node neighbors(i)],[start_node neighbors(i)]); % transition matrix
    Q_power = Q^t; 
    sev_orders(i) = sev0(Q_power);
end
[~,IX]=sort(sev_orders,'descend');
neighbors=neighbors(IX);
%%
for i=1:length(neighbors)
    curr_neighbor = neighbors(i);
    [commrun commsize sev] = sev_node_pass ([start_node curr_neighbor], ...
        adj_matrix, t, max_size);
    %commrun
    if ((ismember(start_node(1),commrun)))% &&(commsize~=graph_size))
        break
    end
end

end

function [ commrun commsize sev ] = sev_node_pass( start_node, adj_matrix, t, max_size)
%SEV_NODE Greedy optimisation for the sev-community of a node, with EM
%   Doesn't necessary give you a community with the original node in it!
%   Copyright (c) 2010 Yun William Yu

graph_size = length(adj_matrix);

%A = adj_matrix(community,community);
D = diag(sum(adj_matrix,2)); % diagonal degree matrix
P = (D^-1)*adj_matrix;  % Markovian transition matrix
%P = full(P);

commrun = zeros(1,graph_size); % Initialise family of communities
commsize = zeros(1,1); % This is really duplicate information to commrun
sev = zeros(1,1); % Vector of severabilities for each community

%disp(t)
community = start_node;
sev_1 = -100000;
%global I2;
j=0;
while length(community)<max_size
    j=j+1;
    if true % (mod(j,1)==0)
        if j~=1
            fprintf('\b\b\b\b')
            %highlighted_image(I2,community);
            %input('press key','s')
        else
            fprintf('.')
        end
        fprintf('%4d',length(community))    % Progress bar
    end
    adj_nodes = find(sum(adj_matrix(community,:),1)); % Only want the adjacent nodes
    neighbours = setdiff(adj_nodes,community);
    % KL-switch every third step
    if (mod(j,3)==0)
        sev_2 = zeros(length(community)+length(neighbours),1);
        % Check the removal of all nodes but the seed node
        parfor i=2:length(community)
            red_comm=community;
            red_comm(i)=[];
            Q=P(red_comm,red_comm);
            sev_2(i)=sev0(Q^t);
        end
        % Check the addition of all possible neighbours
        parfor i=(length(community)+1):length(sev_2)
            Q=P([community neighbours(i-length(community))],[community neighbours(i-length(community))]);
            sev_2(i) = sev0(Q^t);
        end
        % Choose the new community with highest severability
        max_element = find(sev_2==max(sev_2),1,'first');
        if max_element > length(community)
            community=[community neighbours(max_element-length(community))];
        else
            community(max_element)=[];
        end
    else
    % Greedy optimisation on all other steps
        sev_2 = zeros(length(neighbours),1);      
        parfor i=1:length(neighbours)
            Q = P([community neighbours(i)],[community neighbours(i)]);
            sev_2(i) = sev0(Q^t);
        end
        % Choose the new community with highest severability
        max_element = find(sev_2==max(sev_2),1,'first');
        community = [community neighbours(max_element)];
    end

    % If the severability is higher than the current best community, store
    % it; otherwise, we retain the old one.
    if sev_2(max_element)>sev_1
        sev_1=sev_2(max_element);
        comm_final=community;
    end
end
%end

%comm_final
[comm_final sev_1]=sev_steepest_descent(adj_matrix, comm_final, t, max_size);

commrun(1:length(comm_final))=(comm_final);
sev=sev_1;
commsize=length(find(commrun));
    
    

end


