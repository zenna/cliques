function vi=varinfo(partition_vectors,mode)

if nargin < 2
    mode = 'standard';
end

% initialise some constants
number_of_partitions = size(partition_vectors,1);
n = size(partition_vectors,2);
vi_tot=0;
nodes = [1:n];

% The actual VI is the mean of all one to one comparisons..
for i = [1:number_of_partitions]
    % get partition 1 and adjust numbering of nodes
    partition_1 = partition_vectors(i,:);
    partition_1 = double(partition_1)+1;
    % encode community assignments in matrix A_1
    A_1 = sparse(partition_1,nodes,1);
    n_1_all = sum(A_1,2);
    
    % look at all different pairs
    for j = [1:i-1]
        % same procedure as for partition vector 1
        partition_2 = partition_vectors(j,:);
        partition_2 = double(partition_2)+1;
        % order/dimension reversed here
        A_2 = sparse(nodes,partition_2,1);
        n_2_all = sum(A_2,1)';
        
        % n12_all == frequency of nodes occuring in community i,j in
        % partition 1,2
        n_12_all = A_1*A_2;
        
        % only take nonzero values into account...
        [rows,cols,n_12] = find(n_12_all);
        
        n_1 = n_1_all(rows);
        n_2 = n_2_all(cols);
        
        
        % compute variation of information as (H(C|C') + H(C'|C))
        % normalized by log(n)
        vi = sum(n_12.*log(n_12.^2./(n_1.*n_2)));
        if strcmp(mode,'standard')
            vi = -1/(n*log(n))*vi;
            
        elseif strcmp(mode,'renormalized')
            % renormalized by joint entropy, written in a way n and some
            % minus signs cancel out
            H_ClChat = sum(n_12.*log(n_12/n));
            vi = vi/(H_ClChat);
        end
        
        % accumlate
        vi_tot=vi_tot+vi;
        
    end
end

vi = 2*vi_tot/(number_of_partitions*(number_of_partitions-1));

end