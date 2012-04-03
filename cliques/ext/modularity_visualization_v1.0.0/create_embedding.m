function h=create_embedding(filename)

[q_list, partition_vectors] = parse_file(filename);
%Following fails when there is a partition of all nodes in 0
vi_matrix = calculate_vi_matrix(partition_vectors);

%n = 4
%vi_matrix = abs(gallery('randcorr',n).*((gallery('randcorr',n) ~= 1)))
%q_list = rand(n,1)

%Following fails when there is 0 in the vi_matrix (i.e. same cluster
%appearing twice, failure of modularity measure)
h = plotSpace(cca(vi_matrix,2,5000),q_list)

end