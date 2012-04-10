%Test the uniform sampling output. Run the sample_landscape script and 
%load sampling file into sampling_test.

num_nodes = size(sampling_test,2);
numbers = num_nodes.^[0:num_nodes-1];

test = sampling_test*numbers'; % create unique id for each partition

% check if number of partitions found equals bell number
bell_nr6 = 203
nr_partitions_found = length(unique(test))
% get all unique ids
partition_ids = unique(test);

% get distribution by checking all ids
distribution = zeros(1,nr_partitions_found);
for i=1:nr_partitions_found
    distribution(i) = nnz(test == partition_ids(i));
end

% plot normalized distribution
plot(distribution/length(test));

% compare to perfect uniform sampling
sampling_uniform = randi(nr_partitions_found,size(test));
% sampling_uniform = randi(nr_partitions_found,1,1e8);
distribution_uni = zeros(1,nr_partitions_found);
for i=1:nr_partitions_found
    distribution_uni(i) = nnz(sampling_uniform == i);
end
figure
plot(distribution_uni/length(sampling_uniform))