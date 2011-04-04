function [ commrun commsize sev ] = sev_node( start_node, adj_matrix, t, push_ahead, radius)
%SEV_NODE Greedy optimisation for the sev-community of a node, with EM
%   Will always give you a community with the original node in it.
%   Copyright (c) 2010 Yun William Yu
commrun=[];
%commsize=0;
%sev=0;

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
    Q_power = Q^t; %quick_power([community adj_nodes(i)],Q,t);
    sev_orders(i) = sev0(Q_power);
end
[~,IX]=sort(sev_orders,'descend');
neighbors=neighbors(IX);
%%
for i=1:length(neighbors)
    curr_neighbor = neighbors(i);
    [commrun commsize sev] = sev_node_pass ([start_node curr_neighbor], ...
        adj_matrix, t, push_ahead, radius);
    %commrun
    if ((ismember(start_node(1),commrun)))% &&(commsize~=graph_size))
        break
    end
end

end

function [ commrun commsize sev ] = sev_node_pass( start_node, adj_matrix, t, push_ahead, radius)
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
%for j=1:(graph_size-length(start_node))
%global I2;
for j=1:(radius-length(start_node))
    if (mod(j,10)==0)
        if j~=10
            fprintf('\b\b\b\b')
            %highlighted_image(I2,community);
            %input('press key','s')
        else
            fprintf('.')
        end
        fprintf('%4d',j)    % Progress bar
    end
    adj_nodes = find(sum(adj_matrix(community,:),1)); % Only want the adjacent nodes
    sev_2 = zeros(length(adj_nodes),1);      
    parfor i=1:length(adj_nodes)
        if isempty(find(community==adj_nodes(i),1,'first'))
            Q = P([community adj_nodes(i)],[community adj_nodes(i)]); % transition matrix
            Q_power = Q^t; %quick_power([community adj_nodes(i)],Q,t);
            sev_2(i) = sev0(Q_power);
        end
    end

    % Choose the new community with highest severability
    max_element = find(sev_2==max(sev_2),1,'first');
    community = [community adj_nodes(max_element)];

    % If the severability is higher than the current best community, store
    % it; otherwise, we retain the old one.
    %disp(sev_2(max_element))
    if sev_2(max_element)>sev_1
        sev_1=sev_2(max_element);
        comm_final=community;
        push_counter = 0;
    % Not only retain the old one, but also end the process at the local
    % minimum found if push_ahead isn't set.
    else
        push_counter = push_counter + 1;
        if push_counter==push_ahead
            break
        end
    end
end
%end

%comm_final
[comm_final sev_1]=sev_steepest_descent(adj_matrix, comm_final, t, radius);

commrun(1:length(comm_final))=(comm_final);
sev=sev_1;
commsize=length(find(commrun));
    
    

end


