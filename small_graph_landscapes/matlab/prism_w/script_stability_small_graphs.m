% script to run a stability analysis on the "small graphs" for landscape
% analysis.
weights = 1:.01:1.5;
time = logspace(log10(0.5),1,100);

%%
for e=weights
    name = ['prism_w' num2str(e)];
    A = load([name '.edj']);
    
    % create graph from file
    nr_nodes = max(max(A(:,1:2)))+1;
    A = sparse(A(:,1)+1,A(:,2)+1,A(:,3),nr_nodes,nr_nodes);
    A = A+A';
    
    stability_new(A,time,'out',name);
end 

%%

time_c2 = ones(size(weights));
time_c6 = time_c2; time_c3 = time_c2;
i=1;
vi_matrix = zeros(length(weights),length(time));
for e= weights
    name = ['Stability_' 'prism_w' num2str(e) '.mat'];
    load(name);
    
    vi_matrix(i,:)= VI;
    
    test = find(N==6,1,'last');
    if ~(isempty(test))
    time_c6(i) = test;
    end
    
    test = find(N==3,1,'last');
    if ~isempty(test)
    time_c3(i) = test;
    end
    
    test = find(N==2,1,'first');
    if ~isempty(test)
    time_c2(i) = test;
    end
    
    i = i+1;

    %     hold all
    %     script_plot_correlation_stability
end