function[cov_adj]=cov_adjacency(cov_file,top_file,covtable_file)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Stefano Meliga
% Institution: Imperial College London
% Project: Graph Clustering of Atomic Networks for Protein Dynamics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This funcion builds the binding energy matrix from FIRST
% covalent bonds output
% INPUTS: cov_file(FIRST's cov.out), top_file (GROMACS' topology),
% covtable_file
% 
% OUTPUTS:
% 
% INTERNAL VARIABLES:
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MAXDIM=5000;
%inizialization matrix
%A=zeros(MAXDIM,MAXDIM);
cov_matrix=sparse(MAXDIM,MAXDIM);
%loading files
cov_edges=load(cov_file);
covlist=load(covtable_file); % in kJ/mol
covtable=zeros(5);
%length(covlist)
for i=1:length(covlist)
    covtable(covlist(i,1),covlist(i,2))=covlist(i,3);
end
covtable
%LOAD FILE .top
file_id=fopen(top_file);
i=1;
%skip headers and go to atoms section
while isempty(strfind(fgetl(file_id),'[ atoms ]'))==1
end
fgetl(file_id);

for i=1:MAXDIM     %start loading
    top_cell(i) = {fgetl(file_id)};
    if isempty(strfind(top_cell{i},'[ bonds ]'))==0
        N=i-2
        i=MAXDIM;
    end
    i=i+1;
end
fclose(file_id);

for b=1:size(cov_edges,1) % b = bond %CHECK LAST LINE OF THE FILE!!!
    %atom1=cov_edges(b,1);
    %atom2=cov_edges(b,2);
    pair=[0;0];
    oplsdata = ['opls_0';'opls_0'];
    opls = cellstr(oplsdata);
    %opls(1)= {'opls_0'};
    %opls(2)={'opls_0'};
    %b
    for atom=1:2
        % acquisition fields: type and atom
        top_cell{cov_edges(b,atom)};
        [token, rem]=strtok(top_cell{cov_edges(b,atom)});
        [opls{atom} rem]=strtok(rem); % atom type opls_xxx* %cell array
        [token, rem]=strtok(rem);
        [token, rem]=strtok(rem);
        name(atom)={strtok(rem)}; %I want it to be a cell array of strings
        % getting element 
        element(atom)=name{atom}(1); %first character only
        %mapping for cov_table
        switch(element(atom))
            case 'H'
                pair(atom)=1;
            case 'S'
                pair(atom)=2;
            case 'N'
                pair(atom)=3;
            case 'C'
                pair(atom)=4;
            case 'O'
                pair(atom)=5;
            otherwise
                error('element not recognized')
        end
    end
    
    if covtable(pair(1),pair(2)) > 2
        error('Unknown bond type')
    elseif covtable(pair(1),pair(2)) < -2
        cov_matrix(cov_edges(b,1),cov_edges(b,2)) = covtable(pair(1),pair(2));
        cov_matrix(cov_edges(b,2),cov_edges(b,1)) = covtable(pair(1),pair(2));
    else
        %shift=abs(cov_edges(b,2)-cov_edges(b,1)); %distance of atom indexes
        switch abs(covtable(pair(1),pair(2)))
            case {0} %C?C
                if strcmp(opls{1},'opls_145') && (strcmp(opls{2},'opls_145') || strcmp(opls{2},'opls_166'))
                    cov_matrix(cov_edges(b,1),cov_edges(b,2)) = -(602+346)/2; %kJ/mol %1 + 1/2 bond
                    cov_matrix(cov_edges(b,2),cov_edges(b,1)) = -(602+346)/2; %kJ/mol %1 + 1/2 bond                    
                elseif strcmp(opls{1},'opls_507') && strcmp(opls{2},'opls_508')% only His (and maybe Trp), 507 always before 508
                    cov_matrix(cov_edges(b,1),cov_edges(b,2)) = -602; %kJ/mol %double bond
                    cov_matrix(cov_edges(b,2),cov_edges(b,1)) = -602; %kJ/mol %double bond
                else
                    cov_matrix(cov_edges(b,1),cov_edges(b,2)) = -346; %kJ/mol %single bond
                    cov_matrix(cov_edges(b,2),cov_edges(b,1)) = -346; %kJ/mol %single bond
                end
            case {1} %C?N
                if strcmp(opls{1},'opls_302') && strcmp(opls{2},'opls_300') %NOTE: the only case is Arg and C always before N
                    cov_matrix(cov_edges(b,1),cov_edges(b,2)) = -(615+305)/2; %kJ/mol %1 + 1/2 bond
                    cov_matrix(cov_edges(b,2),cov_edges(b,1)) = -(615+305)/2; %kJ/mol %1 + 1/2 bond
                elseif strcmp(opls{1},'opls_511') && strcmp(opls{2},'opls_506')% only His, 511 always before 506
                    cov_matrix(cov_edges(b,1),cov_edges(b,2)) = -615; %kJ/mol %double bond
                    cov_matrix(cov_edges(b,2),cov_edges(b,1)) = -615; %kJ/mol %double bond
                else
                    cov_matrix(cov_edges(b,1),cov_edges(b,2)) = -305; %kJ/mol %single bond
                    cov_matrix(cov_edges(b,2),cov_edges(b,1)) = -305; %kJ/mol %single bond
                end
            case {2} %C?O
                if strcmp(opls{1},'opls_271') && strcmp(opls{2},'opls_272') %only Asp and Glu and C always before O
                    cov_matrix(cov_edges(b,1),cov_edges(b,2)) = -(799+358)/2; %kJ/mol %1 + 1/2 bond
                    cov_matrix(cov_edges(b,2),cov_edges(b,1)) = -(799+358)/2; %kJ/mol %1 + 1/2 bond
                elseif strcmp(opls{1},'opls_235') && strcmp(opls{2},'opls_236')% anywhere, 235 always before 236
                    cov_matrix(cov_edges(b,1),cov_edges(b,2)) = -799; %kJ/mol %double bond
                    cov_matrix(cov_edges(b,2),cov_edges(b,1)) = -799; %kJ/mol %double bond
                else
                    cov_matrix(cov_edges(b,1),cov_edges(b,2)) = -358; %kJ/mol %single bond
                    cov_matrix(cov_edges(b,2),cov_edges(b,1)) = -358; %kJ/mol %single bond
                end
            otherwise ('maths breakdown!')        
        end
    end
end
if length(cov_edges)~= nnz(cov_matrix)/2
    cov_matrix
    length(cov_edges)
    nnz(cov_matrix)/2
    error('COVALENT BONDS: edges number and bonds number do not match')
end
cov_adj=sparse(round(-120*cov_matrix/6.022)); %J/m^2 (SI unit)
%check isolated sites
M=max(cov_adj(1:N,1:N));
j=0;
if min(M)==0
    warning('There are isolated sites:')
    for i=1:N
        if M(i)==0
            j=j+1;
            isol(j)=i;
        end
    end
    isol;
end
