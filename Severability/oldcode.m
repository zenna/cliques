%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sev0.m 2010-11-23
        % Calculate the distance (abs) between evolution of the
        % uniform distribution to the quasi-stationary one
        %[V ~]=eigs(Q',2);
        %q_dist=V(:,1)'/norm(V(:,1));  % Finds the quasi-stationary distribution of Q (normalised)
        %unif_evol=sum(Q_power,1); % Evolution of the uniform distribution by time
        %unif_evol=unif_evol/norm(unif_evol); % Normalised
        
        %start_dist = ones(1,community_size)/norm(ones(1,community_size));        
        %mixing_distance = 1-abs(q_dist*unif_evol')/(norm(q_dist)*norm(unif_evol));
        %mixing_benchmark = abs(q_dist*start_dist')/(norm(q_dist)*norm(start_dist));
        %distance = min([pdist([q_dist; unif_evol]) pdist([-q_dist; unif_evol])]);
        %start_dist = ones(1,community_size)/norm(ones(1,community_size));
        %distance_benchmark = min([pdist([q_dist; start_dist]) pdist([-q_dist; start_dist])]);
        %%mixing = log(distance)-log(distance_benchmark);
        %real_mix=log(mixing_distance);

        %mixing = 1-mean(pdist(Q_power,'cosine'));
        
        % The above line is too time-consuming, but by Barhum, Goldreich,
        % and Shraibman, 2007, we can approximate the average distance by
        % taking a random sampling of the combinatorial possibilities.
%         cos_dist_vec=zeros(1,10*community_size);
%         for i=1:length(cos_dist_vec)
%             cos_dist_vec(i)=pdist(Q_power(randi(community_size,1,2),:),'cosine');
%         end
%        mixing = 1-mean(cos_dist_vec);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sev_node.m 2010-11-24
% Hash table (ugly hack) of the communities and their matrix powers
% key_table={};
% key_time_table=[];
% power_table={};
%     function [Q_quick_power] = quick_power(comm,Q,t)
%         comm_key = keygen(comm);
%         location=find(strcmp(key_table,comm_key),1,'first');
%         if isempty(location)
%             location=length(key_table)+1;
%             key_table{location}=comm_key;
%             Q_quick_power = Q^t;
%             power_table{location}=Q_quick_power;
%             key_time_table(location)=t;
%         else
%             %Q_quick_power = power_table{location}*Q^(t-key_time_table(location));
%             Q_quick_power = Q^t;
%             power_table{location}=Q_quick_power;
% %             if Q_quick_power~=Q^t
% %                 disp('error!!!');
% %                 t
% %                 i
% %                 j
% %             end
%             key_time_table(location)=t;
%         end
%     end
%    function [Q_quick_power] = quick_power(comm,Q,t)
%        Q_quick_power=Q^t;
%    end

    %     % Choose the new community with highest severability, if it's greater
    %     % than the previous community's severability, OR if the previous
    %     % community had a severability beneath the threshhold level.
    %     if (sev_2 > sev) || sev < 0.01
    %         sev=sev_2;
    %         r3 = r2;
    %         m3 = m2;
    %         community = [community new_node];
    %     % Else, end the loop.
    %     else
    %         disp('Local severability maximum found')
    %         break
    %     end




