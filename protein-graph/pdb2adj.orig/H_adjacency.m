function[hb_adj]=H_adjacency(hb_file, sb_file)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Stefano Meliga
% Institution: Imperial College London
% Project: Graph Clustering of Atomic Networks for Protein Dynamics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This funcion builds the weighted adjacency matrix from FIRST
% hydrogen bonds output
% INPUTS: hb_file(FIRST's hbonds.out), sb_file(list of salt bridges obtained processing
% FIRST's *_results.txt)
% 
% OUTPUTS:
% 
% INTERNAL VARIABLES:
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MAXDIM=5000;
Vs=8 %kcal/mol
Rs=3.2 %A
x=0.375 %A

%inizialization matrix
%A=zeros(MAXDIM,MAXDIM);
hb_matrix=sparse(MAXDIM,MAXDIM);
%loading files
hb_list=load(hb_file);
hb=0;
for i=1:length(hb_list)
    if hb_list(i,3)<-0.01 % threshold to ensure integer weigths: note k=-120*E*4.184/6.022
        hb_matrix(hb_list(i,1),hb_list(i,2))=hb_list(i,3); % Kcal/mol
        hb_matrix(hb_list(i,2),hb_list(i,1))=hb_list(i,3);
        hb=hb+1;
    end
end

% matrix of uncharged hydrogen bonds force constants
hb_adj=round(-120*hb_matrix*4.184/6.022); %conversion from kcal/mol to J/A^2

sb_list=load(sb_file);
% calculation of unique force constant for salt bridges
k_sb=round(120*Vs*((Rs+x)/Rs)*(4.184/6.022))
for i=1:length(sb_list)
    i;
    hb_adj(sb_list(i,1),sb_list(i,2))=k_sb;
    hb_adj(sb_list(i,2),sb_list(i,1))=k_sb;
end
if hb ~= nnz(hb_adj)/2
    %hb_adj
        hb
        nnz(hb_adj)/2
    warning('HYDROGEN BONDS: edges number and bonds number do not match')
end
