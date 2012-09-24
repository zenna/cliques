% script to run a stability analysis on the "small graphs" for landscape
% analysis.

time = logspace(-1,1,200);

for e=.1:.2:4
    name = ['4node_w' num2str(e)];
    A = load([name '.edj']);
    
    % create graph from file
    nr_nodes = max(max(A(:,1:2)))+1;
    A = sparse(A(:,1)+1,A(:,2)+1,A(:,3),nr_nodes,nr_nodes);
    A = A+A';
    
    stability_new(A,time,'out',name);
end 