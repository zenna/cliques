function [ orphans ] = findorphans( prospectives, A, t, max_size )
%FINDORPHANS Checks to see if a set of nodes are orphans
%   [ orphans ] = findorphans(A, prospectives)
%   
%   prospectives: [ n1, n2, ...], where {n} are the nodes to be questioned
%              A: adjacency matrix
%              t: Markov time
%       max_size: Limit for the search algorithm
%
%        orphans: all the orphaned nodes (as a vector)
%
%   Copyright (c) 2011 Yun William YU

if any(prospectives<1)||any(prospectives>length(A))
    error('Invalid nodes in "prospectives".');
end

orphans=[];
for i=1:length(prospectives)
    [commrun commsize sev]=sev_node(prospectives(i), A, t, max_size);
    if any(commrun==prospectives(i)) && commsize>1
        orphans(end+1)=prospectives(i);
    end
end



end

