function[hphobes_adj]=hphobes_adjacency(hphobes5A_file, hphobes8A_file)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Stefano Meliga
% Institution: Imperial College London
% Project: Graph Clustering of Atomic Networks for Protein Dynamics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This funcion builds the weighted adjacency matrix from two FIRST outputs:
% the lists of hydrophobic bonds with thresholds 5 and 8 Angstrom.
% The well depths are taken from the article Hydrophobic Potential of Mean
% Force as... by M.S. Lin, June 2007, Cell Press.
% INPUTS: hphobes5A_file(FIRST's hphobes.out@5A), hphobes8A_file(FIRST's
% hphobes.out@8A)
% 
% OUTPUTS:
% 
% INTERNAL VARIABLES:
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MAXDIM=5000;

%inizialization matrix
%A=zeros(MAXDIM,MAXDIM);
hphobes_matrix=sparse(MAXDIM,MAXDIM);

%loading files
hphobes5A_list=load(hphobes5A_file);
hphobes8A_list=load(hphobes8A_file);
%{
for i=1:length(hphobes8A_list)
    hphobes_matrix(hphobes8A_list(i,1),hphobes8A_list(i,2))=-0.1; % kcal/mol
    hphobes_matrix(hphobes8A_list(i,2),hphobes8A_list(i,1))=-0.1; % kcal/mol
end
for i=1:length(hphobes5A_list)
    if hphobes_matrix(hphobes5A_list(i,1),hphobes5A_list(i,2))==0
        warning()
    hphobes_matrix(hphobes5A_list(i,1),hphobes5A_list(i,2))=hphobes_matrix(hphobes5A_list(i,1),hphobes5A_list(i,2))-0.7; % kcal/mol
    hphobes_matrix(hphobes5A_list(i,2),hphobes5A_list(i,1))=hphobes_matrix(hphobes5A_list(i,2),hphobes5A_list(i,1))-0.7; % kcal/mol
end
%}

for i=1:length(hphobes5A_list)
    hphobes_matrix(hphobes5A_list(i,1),hphobes5A_list(i,2))=-0.8; % kcal/mol
    hphobes_matrix(hphobes5A_list(i,2),hphobes5A_list(i,1))=-0.8; % kcal/mol
end
c=0;
for i=1:length(hphobes8A_list)
    if hphobes_matrix(hphobes8A_list(i,1),hphobes8A_list(i,2))==0
        hphobes_matrix(hphobes8A_list(i,1),hphobes8A_list(i,2))=-0.1; % kcal/mol
        hphobes_matrix(hphobes8A_list(i,2),hphobes8A_list(i,1))=-0.1; % kcal/mol
        c=c+1;
    end
end

% matrix of uncharged hydrophobic bonds force constants
hphobes_adj=round(-120*hphobes_matrix*4.184/6.022); %conversion from kcal/mol to J/A^2

if length(hphobes5A_list)+c ~= nnz(hphobes_adj)/2
    hphobes_adj
    length(hphobes5A_list)+c
    nnz(hphobes_adj)/2
    error('HYDROPHOBIC BONDS: edges number and bonds number do not match')
end
