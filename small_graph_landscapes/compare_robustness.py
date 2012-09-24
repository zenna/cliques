import getopt, sys
import simplejson
from collections import defaultdict
import matplotlib.pyplot as plt
import math
import IPython

# TODO

def usage():
    print """-x input prefix e.g. ~/path_graph (where files ~/path_graph_coords.mat etc have been generated
    or you can be specific (will override prefix):
    -e edges_file
    -r energy_file (e.g. stability)
    -c coordinates
    -g graph_file
    -p partitions
    -o output_file
    -b basin_file
    """

def parse_args():
    """Parse arguments of script"""
    try:
        opts, args = getopt.getopt(sys.argv[1:], "b:g:p:r:c:e:h:x:", ["help", "output="])
    except getopt.GetoptError, err:
        print str(err)
        usage()
        sys.exit(2)
    
    # default arguments
    coords_file = None
    edges_file = None
    energy_file = None
    partitions_file = None
    graph_file = None
    basins_file = None
    input_file = None
    
    # assign file names from standard prefix
    for o, a in opts:
        if o in ("-x", "--prefix"):
            coords_file = '%s_%s.mat'%(a, 'coords')
            edges_file = '%s_%s.edj'%(a, 'landscape_edgelist')
            energy_file = '%s_%s.mat'%(a, 'energy')
            basins_file = '%s_%s.mat'%(a, 'greedy_basins')
            partitions_file = '%s_%s.mat'%(a, 'partitions')
            graph_file = '%s_%s.edj'%(a, 'graph_edgelist')
    
    # in case specific files are given use these
    for o, a in opts:
        if o in ("-c", "--coordinates"):
            coords_file = a
        elif o in ("-e", "--edges"):
            edges_file = a
        elif o in ("-r", "--energy"):
            energy_file = a
        elif o in ("-b", "--basins"):
            basins_file = a
        elif o in ("-p", "--partitions"):
            partitions_file = a
        elif o in ("-g", "--graph"):
            graph_file = a
        elif o in ("-h", "--help"):
            usage()
            sys.exit()
   
    return coords_file, edges_file, energy_file, partitions_file, graph_file, basins_file

def file_to_nested_list(filename, cast_type):
    """Convert file to nested list"""
    f = open(filename, 'r')
    nested_list = []
    for line in f:
        if line:
            l = [cast_type(x) for x in line.split() if x]
            nested_list.append(l)
            
    f.close()
    return nested_list

def file_to_energy_process(filename, cast_type):
    """Convert energy file to list of energy values. Format: energies[time_index]['values'/'time']"""
    f = open(filename, 'r')
    data = []
    for line in f:
        if line:
            l = [cast_type(x) for x in line.split() if x]
            data_point = {'time':l[0], 'values':l[1:]}
            data.append(data_point)
            
    f.close()
    return data

def file_to_basin_sizes(filename, energies, time_limit):
    """return basin list with basins[basin_id]['time/psize/size']"""
    f = open(filename, 'r')
    data = []
    # a dictionary with default value dictionary with default value list
    basin_to_values = defaultdict(lambda: defaultdict(list))
    time_index = 0
    for line in f:
        if line:
            l = [x for x in line.strip().split(" ") if x]
            time = float(l[0])
            if  time < time_limit:
                # find corresponing energy entries
                if time != energies[time_index]['time']:
                    time_index += 1
                    
                basin_size = 0
                basin_prob_size = 0.0
                # first part is the node, second part the prob to go to basin
                for index, val in enumerate(l[2:]):
                    if index % 2 == 0:
                        nodes_id = int(val)
                        basin_size += 1
                    else:
                        energy = energies[time_index]['values'][nodes_id]
                        node_to_basin_prob = float(val)
                        basin_prob_size += node_to_basin_prob #* energy
                basin_id = int(l[1])    
                basin_to_values[basin_id]['time'].append(time)
                basin_to_values[basin_id]['psize'].append(basin_prob_size)
                basin_to_values[basin_id]['size'].append(basin_size)
    
    f.close()
    return basin_to_values

def extract_basins_per_time(filename, time_limit):
    basins_per_time = [[]]
    with open(filename,'r') as f:
        for line in f.readlines():
            if line:
                l = [x for x in line.strip().split(" ") if x]
                time = float(l[0])
                if  time < time_limit:
                    #todo
                        


def main():
    #IPython.embed()
    # parse inputs
    coords_file, edges_file, energy_file, partitions_file, graph_file, basins_file = parse_args()
    time_limit = 200
    if coords_file:
        print "processing coords"
        try:
            coordinate_list = file_to_nested_list(coords_file, float)
        except:
            pass
            
    if edges_file:
        edge_list = file_to_nested_list(edges_file, int)
    
    if energy_file:
        print "processing energy"
        energy_list = file_to_energy_process(energy_file, float)

    if partitions_file:
        try:
            partition_list = file_to_nested_list(partitions_file, int)
        except:
            print "couldnt open partitions"
 
    if graph_file:
        import networkx as nx
        graph = nx.read_edgelist(graph_file, nodetype=int, data=(('weight',float),))
        pos = nx.spring_layout(graph,iterations=500)
        graph_list = {}
        graph_list['coords'] = [x.tolist() for x in pos.values()]
        graph_list['edges'] = file_to_nested_list(graph_file, float)
        #plt.figure()
        #nx.draw(graph)
 
    if basins_file:
        num_partitions = float(len(partition_list))
        print "number of partitions: " +  str(num_partitions)
                    
        basins = file_to_basin_sizes(basins_file, energy_list,time_limit)        
        basins_per_time = extract_basins_per_time(basins_file,time_limit)           
        IPython.embed()

if __name__ == "__main__":
    main()
