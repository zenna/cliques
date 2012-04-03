%script that converts file edgesJoao2.mat given by Joao to matrix (=adjacency of the graph)
%The adjacency matrix is Q.

 %edges contains the edges
 
 %LARGE FILE!!

function [Q,nodes,links]=ConvertFile3W(adj_file,N)
 edgesW=load(adj_file)
 %open edgesJoao2.mat
 %edges=ans.edges

Q=sparse(zeros(N,N));
links=length(edgesW)
for i=1:links
   Q(edgesW(i,1)+1,edgesW(i,2)+1)=edgesW(i,3);
end
nodes=edgesW(i,1)+1
if Q ~= Q'
   error('The adjacency matrix is not symmetrical')
end