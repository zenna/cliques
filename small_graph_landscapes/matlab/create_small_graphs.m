%script to create small graphs for landscapes..


for e=.001:.002:1
%barbell graph with varying strength of "bridge"
A = [0 1 1 0 0 0
     1 0 1 0 0 0
     1 1 0 e 0 0
     0 0 e 0 1 1 
     0 0 0 1 0 1

     0 0 0 1 1 0];
[i j v] = find(tril(A));
dlmwrite(['barbell_bridge' num2str(e) '.edj'],[i-1 j-1 v],'delimiter', ' ')
end

%%
for e=1.0:.01:1.5
%prism graph with varying strength of "bridges"
A = [0 1.3 0.8 0 e 0
     1.3 0 1 0 0 e
     0.8 1 0 e 0 0
     0 0 e 0 0.9 1 
     e 0 0 0.9 0 1
     0 e 0 1 1 0];
[i j v] = find(tril(A));
dlmwrite(['prism_w' num2str(e) '.edj'],[i-1 j-1 v],'delimiter',' ')
end

%%
for e=.8:.025:1.3
%4 node graph with varying strength of "bridges"
A = [0 1 .1 e
     1 0 e .1
     .1 e 0 1
     e .1 1 0];
[i j v] = find(tril(A));
dlmwrite(['4node_w' num2str(e) '.edj'],[i-1 j-1 v],'delimiter',' ')
end
