function convertAdjMatrixToPajekGraph(A,filename)
% Converts adjacency matrix of undirected network into pajek graph of an
% undirected network. Assumes network to be fully connected.
% Inputs:
%           A:          Adjacency matrix of undirected graph
%       
%           filename:   string with filename of output file,
%                       e.g. 'pajekgraph.net'. If file already exists
%                       contents are overwritten
%
% last revision: 25/3/2010 by Michael  

nr_nodes = length(A);
% open file and write name list
fid = fopen(filename,'w+','n','ascii');
fprintf(fid,'*Vertices %i \n',nr_nodes);
for i= 1:nr_nodes
    fprintf(fid,'%i \"%i\" \n', int32(i), int32(i));
end


% print edges
fprintf(fid,'*Edges\n');
for i=1:nr_nodes
    for j=i:nr_nodes
        link = A(i,j);    
        if(link)
            fprintf(fid,'%i %i %e \n', i, j, link);
        end
    end
end

fclose('all');

end