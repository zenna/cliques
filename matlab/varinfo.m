function [vi,vi_mat] = varinfo(partition_vectors,ComputeParallel)

number_of_partitions = size(partition_vectors,1);    
n = size(partition_vectors,2);
vi_mat = zeros(number_of_partitions);
vi=0;
    
% If all the partitions are identical, vi=0 and there is no need to do the
% rest of the calculations which are computationally expensive.
if  all(all(partition_vectors==repmat(partition_vectors(1,:),number_of_partitions,1)))
    return;
% else
%     % Otherwise, reorder numbering of the communities such that, ie 
%     % [0,0,0,1,1,1] and [1,1,1,0,0,0] become identical.
%     for i=1:number_of_partitions
%         [a,b] = unique(partition_vectors(i,:),'First');
%         [I,I] = sort(b);
%         a=a(I);
%         tmp=ones(1,length(partition_vectors(i,:)))*(max(partition_vectors(i,:))+1);
%         for j=1:length(a)
%             tmp(partition_vectors(i,:)==(a(j)))=j-1;
%         end
%         partition_vectors(i,:)=tmp;
%     end
end

% Select only the partitions which are different 
[partition_vectors,b,c] = unique(partition_vectors,'rows');

number_of_partitions=length(b);

vi_mat = zeros(number_of_partitions);

vi_tot=0;
nodes = 1:n;

if nargin==2 && ComputeParallel
    parfor i = 1:number_of_partitions
        partition_1 = partition_vectors(i,:);
        partition_1 = double(partition_1)+1;
        A_1 = sparse(partition_1,nodes,1);
        n_1_all = sum(A_1,2);
        vi_mat_row=vi_mat(i,:);

        for j = 1:i-1
            partition_2 = partition_vectors(j,:);
            partition_2 = double(partition_2)+1;
            A_2 = sparse(nodes,partition_2,1);
            n_2_all = sum(A_2,1)';
            n_12_all = A_1*A_2;

            [rows,cols,n_12] = find(n_12_all);

            n_1 = n_1_all(rows);
            n_2 = n_2_all(cols);

            vi = sum(n_12.*log(n_12.^2./(n_1.*n_2)));
            vi = -1/(n*log(n))*vi;
            
            vi_mat_row(j)=vi;
            
            vi_tot=vi_tot+vi;

        end
        vi_mat(i,:)=vi_mat_row;
    end
else
    for i = 1:number_of_partitions
        partition_1 = partition_vectors(i,:);
        partition_1 = double(partition_1)+1;
        A_1 = sparse(partition_1,nodes,1);
        n_1_all = sum(A_1,2);

        for j = 1:i-1
            partition_2 = partition_vectors(j,:);
            partition_2 = double(partition_2)+1;
            A_2 = sparse(nodes,partition_2,1);
            n_2_all = sum(A_2,1)';
            n_12_all = A_1*A_2;

            [rows,cols,n_12] = find(n_12_all);

            n_1 = n_1_all(rows);
            n_2 = n_2_all(cols);

            vi = sum(n_12.*log(n_12.^2./(n_1.*n_2)));
            vi = -1/(n*log(n))*vi;
            vi_mat(i,j)=vi;
            vi_tot=vi_tot+vi;

        end
    end
end

vi_mat_new = zeros(number_of_partitions,length(c));

for i=1:number_of_partitions
    vi_mat_new(i,:) = vi_mat(i,c);
end
vi_mat_new=vi_mat_new(c,:);

vi_mat = vi_mat_new+vi_mat_new';%vi_mat+vi_mat';

vi = mean(squareform(vi_mat));

end
