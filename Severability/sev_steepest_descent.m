function [community sev] = sev_steepest_descent(A, community, t, max_size)
%SEV_STEEPEST_DESCENT Finds the local maximum of severability
%   A = adjacency matrix
%   community = [v_1 v_2 ... v_k], where v_1 is the number of the node
%   t = scalar integer time
%   sev = scalar integer


% Doesn't make sense for the max size to be greater than the community size
max_size = min(max_size,size(A,1));

P=diag(sum(A,2).^-1)*A;
sev_1 = sev0(P(community,community)^t);
increase=true;
firsttime=true;
while increase==true
    if true
        fprintf('\b\b\b\b')
        fprintf('%4d',length(community))    % Progress bar
    end
    adj_nodes=find(sum(A(community,:),1)); % Get all adjacent nodes
    length_adj_nodes=length(adj_nodes);
    sev_2 = zeros(length_adj_nodes+length(community),1);
    parfor i=1:length_adj_nodes
        if isempty(find(community==adj_nodes(i),1,'first'))
            Q=P([community adj_nodes(i)],[community adj_nodes(i)]);
            Q_power = Q^t;
            sev_2(i) = sev0(Q_power);
        end
    end
    parfor i=1:length(community)
        % Reduce the community
        red_comm=community;
        red_comm(i)=[];
        Q=P(red_comm,red_comm);
        Q_power = Q^t;
        sev_2(i+length_adj_nodes) = sev0(Q_power);
    end
    
    % Choose the new community with highest severability
    max_element = find(sev_2==max(sev_2),1,'first');
    if sev_2(max_element) > sev_1
        increase = true;
        sev_1=sev_2(max_element);
        if max_element > length(adj_nodes)
            community(max_element - length(adj_nodes))=[];
        else
            community=[community adj_nodes(max_element)];
        end
    else
        increase = false;
    end
    % Stops the process if the community goes beyond the max_size
    if (length(community)>max_size)
        increase = false;
        %community =community(1);
    elseif (length(community)==1)
        increase = false;
    end
end

sev = sev_1;

end



