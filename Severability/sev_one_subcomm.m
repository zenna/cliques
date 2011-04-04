function [ subcomm sev ] = sev_one_subcomm( start_node, adj_matrix, t )
%SEV_NODE Greedy optimisation for the sev-community of a node
%   Stops when the severability no longer increases, so only one
%   "subcommunity"
%   Copyright (c) 2010 Yun William Yu

graph_size = length(adj_matrix);

%A = adj_matrix(community,community);
D = diag(sum(adj_matrix),2); % diagonal degree matrix
P = (D^-1)*adj_matrix;  % Markovian transition matrix
%P = full(P);

t;
disp(t)
community = start_node;
sev_1 = -100000;
for j=1:graph_size
    if (mod(j,10)==0)
        fprintf('\b\b\b\b')
        fprintf('%4d',j)    % Progress bar
    end
    adj_nodes = find(sum(adj_matrix(community,:),1)); % Only want the adjacent nodes
    sev_2 = zeros(length(adj_nodes),1);      
    for i=1:length(adj_nodes)
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
    if sev_2(max_element)>sev_1
        sev_1=sev_2(max_element);
        comm_final=community;
    % Not only retain the old one, but also end the process at the local
    % minimum found if push_ahead isn't set.
    else
        break
    end



end
subcomm=comm_final;
sev=sev_1;
   
    

end


