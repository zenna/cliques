function[adj]=adjacency(cov_file,top_file,covtable_file,hb_file,sb_file,hphobes5A_file,hphobes8A_file)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Stefano Meliga
% Institution: Imperial College London
% Project: Graph Clustering of Atomic Networks for Protein Dynamics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This funcion builds the final weighted adjacency matrix using the functions:
% cov_adjacency, H_adjacency, hphobes_adjacency
% 
% INPUTS: just the inputs needed by the internal functions
% 
% OUTPUTS:
% 
% INTERNAL VARIABLES:
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

adj = cov_adjacency(cov_file,top_file,covtable_file) + H_adjacency(hb_file, sb_file) + hphobes_adjacency(hphobes5A_file, hphobes8A_file);

if adj ~= adj'
    error('The adjacency matrix is not symmetric')
else
    disp('Total number of edges:' )
    edges = nnz(adj)/2

    %%plots
    %spy(adj)
    %figure; imagesc (adj(1:100,1:100)); figure(gcf)
    %figure; imagesc (adj(700:900,2800:3000)); figure(gcf)

    %colormapeditor
    %select figure, modify scales, click on colorbar on figure, save images
    w_matrix=zeros(edges,3);
    e=0;
    for i=1:length(adj)
        for j=1:length(adj)
            if adj(i,j)~=0
                e=e+1;
                w_matrix(e,1)=i-1; %-1 to avoid shift of node reordering of Louvain's convert function
                w_matrix(e,2)=j-1;
                w_matrix(e,3)=round(adj(i,j)/2+0.1); %half of weigth because it is double by convert function
            end
        end
    end
    dlmwrite('adj_file.txt', w_matrix, '\t');
end

