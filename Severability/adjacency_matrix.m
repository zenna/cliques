function [ A ] = adjacency_matrix( adj_file )
%Adjacency Matrix
%   Generates the adjacency matrix from an adjacency file (either weighted
%   or unweighted).
%   Copyright (c) 2010 Yun William Yu
X = load(adj_file);

%offset = 0;
if (X(1,1)==0)
    offset = 1;
else
    offset = 0;
end

% Generates the adjacency matrix A
num_vertices = max([X(:,1);X(:,2)]) + offset;
A = sparse(num_vertices,num_vertices);
if (length(X(1,:))==3)
    for i = 1:length(X(:,1))
        A(X(i,1)+offset, X(i,2)+offset) = A(X(i,1)+offset, X(i,2)+offset)  + X(i,3);
    end
else
    for i = 1:length(X(:,1))
        A(X(i,1)+offset, X(i,2)+offset) = A(X(i,1)+offset, X(i,2)+offset) + 1;
    end
end

end

