%script to create landscapes graphs, project them and get the basins.

DIR_REP = '/home/mts09/repositories';
DIR_BIN = '/home/mts09/PhD/ext_programmes_and_lib/cliques/';

weights = 1:.01:1.5;


for e=weights
    name = ['prism_w' num2str(e)];
    path_to_graph = [pwd() name];
    unix(sprintf('%s/tests/test_louvain %s',DIR_BIN, path_to_graph));
end