import getopt, sys
from collections import defaultdict
import matplotlib.pyplot as plt
from matplotlib import cm
import math
import IPython
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

# TODO

def parse_args():
    """Parse arguments of script"""
    try:
        opts, args = getopt.getopt(sys.argv[1:], "x:")
    except getopt.GetoptError, err:
        print str(err)
        sys.exit(2)
    
    # assign standard prefix
    for o, a in opts:
        if o in ("-x", "--prefix"):
            prefix = a            

    return prefix

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
    # a dictionary with default value: dictionary with default value list
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
    # basins per file is a dictionary of lists dict['time'] = [bid1, bid2..]
    basins_per_time = defaultdict(list)
    time = 0
    with open(filename,'r') as f:
        for line in f.readlines():
            if line:
                l = [x for x in line.strip().split(" ") if x]
                time = float(l[0])
                if time > time_limit:
                    break
                else:
                    # append number of basin to list
                    basins_per_time[time].append(int(l[1]))
    return basins_per_time
                        


def main():
    #IPython.embed()
    # parse inputs, get prefix
    global_prefix = parse_args()
    hierarchy = "_0"
    prefix = "out" + hierarchy
    case = 0    # case 0 == prism_w graph
    if case == 0:
        eps_list = np.arange(1,1.51,0.01).tolist()
    
    # prepare stuff to plot
    nr_basins_time = [[]]
    
    for eps in eps_list:
        graphname = global_prefix.split("/")[1]
        dir = "output_" + graphname + str(eps).strip(".00") + "/"
        coords_file = '%s%s%s_%s.mat'%(global_prefix,dir,prefix, 'coords')
        edges_file = '%s%s%s_%s.edj'%(global_prefix, dir, prefix, 'landscape_edgelist')
        energy_file = '%s%s%s_%s.mat'%(global_prefix, dir, prefix, 'energy')
        basins_file = '%s%s%s_%s.bsn'%(global_prefix, dir, prefix, 'greedy_basins')
        partitions_file = '%s%s%s_%s.mat'%(global_prefix, dir, prefix, 'partitions')
        graph_file = '%s%s%s_%s.edj'%(global_prefix, dir, prefix, 'graph_edgelist')


        time_limit = 200

        coordinate_list = file_to_nested_list(coords_file, float)
        edge_list = file_to_nested_list(edges_file, int)
        energy_list = file_to_energy_process(energy_file, float)
        partition_list = file_to_nested_list(partitions_file, int)
            
        # get graph for analysis -- atm not used 
        #import networkx as nx
        #graph = nx.read_edgelist(graph_file, nodetype=int, data=(('weight',float),))
        #pos = nx.spring_layout(graph,iterations=500)
        #graph_list = {}
        #graph_list['coords'] = [x.tolist() for x in pos.values()]
        #graph_list['edges'] = file_to_nested_list(graph_file, float)
        #plt.figure()
        #nx.draw(graph)


        #####################
        # basin analysis
        num_partitions = float(len(partition_list))
        print "number of partitions: " +  str(num_partitions)
        
        # parse files in two ways
        basins = file_to_basin_sizes(basins_file, energy_list,time_limit)        
        basins_per_time = extract_basins_per_time(basins_file,time_limit) 
        
        # get some simpler statistics
        nr_basins_per_time = []
        time_list = []
        for time in sorted(basins_per_time.keys()):
            time_list.append(time)          
            nr_basins_per_time.append( int(len(basins_per_time[time])) )

        nr_basins_time.append(nr_basins_per_time)  

    
    ############################
    # PLOTS

    # remove empty entry -- why does this exist???    
    del nr_basins_time[0]      
    X, Y  = np.meshgrid(eps_list, time_list)
    fig = plt.figure()
    ax = Axes3D(fig)
    Z = np.array(nr_basins_time).transpose()
    ax.plot_surface(X, Y, Z, cmap=cm.jet)
    plt.show()

    #IPython.embed()

if __name__ == "__main__":
    main()
