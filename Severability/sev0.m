% Computes a sev0, which takes the already exponentiated Q-transition
% submatrix and spits back the severability, retention, and mixing
function [sev retention mixing ] = sev0(Q_power)
    community_size=length(Q_power);
    retention = full((1/community_size)*sum(sum(Q_power)));   % Normalised sum of entries

    %cos_sim_vec=zeros(1,10*community_size);
    % Normalise Q_power so we don't have to do it below
%        nQ_power = full(spdiags(sqrt(sum(Q_power.^2,2)).^-1,0,community_size,community_size)*Q_power);
    %for i=1:length(cos_sim_vec)
    %    s = randi(community_size);
    %    t = randi(community_size);
    %    %cos_sim_vec(i)=(Q_power(s,:)*Q_power(t,:)')/(norm(Q_power(s,:))*norm(Q_power(t,:)));
    %    cos_sim_vec(i)=(nQ_power(s,:)*nQ_power(t,:)');
    %end
    %mixing = mean(cos_sim_vec);
    
%        % Try using a variance like approximation
%        avg_Q=mean(nQ_power);
%        navg_Q=avg_Q/norm(avg_Q);
%        cos_avg=zeros(1,community_size);
%        for i=1:community_size
%            cos_avg(i)=(nQ_power(i,:)*navg_Q');
%        end
%        %variance = mean((cos_avg).^2);

%        cos_dist=acos(cos_avg); %/(pi/2);
%        mixing=1-(mean(cos_dist)/acos(community_size^-0.5)); % Mean geodesic distance
%        mixing=1-mean(cos_dist)/(pi/2); % Mean geodesic distance to geographic centroid
%        mixing2=1-mean(pdist(Q_power,'cosine')); % Original pairwise distance

    % New attempt at defining mixing based on total variation distance (1/2 the 1-norm)
    nQ_power = (spdiags(sum(Q_power,2).^-1,0,community_size,community_size))*Q_power;
    avg_Q = mean(nQ_power); % No need to renormalise because ||avg_Q||_1 = 1 (as a prob dist should)
    % Use the first left eigenvector
    %try
    %    [V ~]=eigs(Q_power',2);
    %    q_dist=V(:,1)'/sum(V(:,1));  % Finds the quasi-stationary distribution of Q (normalised)
    %catch exception
    %    q_dist=avg_Q; % Use the average of the matrix if the eigenvector computation fails
    %    disp('Eigenvalue Error');
    %end
    %pdist2(avg_Q, q_dist,'cityblock')*0.5 % Measures the distance between the average and the quasi-stationary distribution
    var_dist = pdist2(nQ_power, avg_Q,'cityblock')*0.5; % Total variation distance
    variance = (mean((var_dist)));
    mixing = 1-variance;

    % Penalise completely empty matrices severely
    if nnz(nQ_power)==0
        mixing=0;
    elseif ~isempty(find(sum(nQ_power,2)==0,1,'first')) % If at least one row is equal to 0
        mixing=0;
    end

    %sev = retention + mixing;
    sev = (retention+mixing)/2;
    %sev2 = retention + 1 - mean(pdist(nQ_power,'cityblock'))*0.5;
end
