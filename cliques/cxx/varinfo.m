function vi=varinfo(partition_vectors)

number_of_partitions = size(partition_vectors,1);    
n = size(partition_vectors,2);
vi_tot=0;
nodes = [1:n];

for i = [1:number_of_partitions]
    partition_1 = partition_vectors(i,:);
    partition_1 = double(partition_1)+1;
    A_1 = sparse(partition_1,nodes,1);
    n_1_all = sum(A_1,2);
            
    for j = [1:i-1]
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
        
        vi_tot=vi_tot+vi;
                   
    end
end

vi = 2*vi_tot/(number_of_partitions*(number_of_partitions-1));

end