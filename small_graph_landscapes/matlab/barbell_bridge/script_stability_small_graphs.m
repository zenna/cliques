% script to run a stability analysis on the "small graphs" for landscape
% analysis.
clc;
clear

time = logspace(-1,1,200);
weights = .001:.002:1;

for e=weights
    name = ['barbell_bridge' num2str(e)];
    A = load([name '.edj']);
    
    % create graph from file
    nr_nodes = max(max(A(:,1:2)))+1;
    A = sparse(A(:,1)+1,A(:,2)+1,A(:,3),nr_nodes,nr_nodes);
    A = A+A';
    
    stability_new(A,time,'out',name);
end 


%%

time_c2 = ones(size(weights));
time_c6 = time_c2; time_c4 = time_c2;
i=1;
for e= weights
    name = ['Stability_' 'barbell_bridge' num2str(e) '.mat'];
    load(name);
    
    test = find(N==6,1,'last');
    if ~(isempty(test))
    time_c6(i) = test;
    end
    
    test = find(N==4,1,'last');
    if ~isempty(test)
    time_c4(i) = test;
    end
    
    test = find(N==2,1,'first');
    if ~isempty(test)
    time_c2(i) = test;
    end
    
    i = i+1;

    %     hold all
    %     script_plot_correlation_stability
end