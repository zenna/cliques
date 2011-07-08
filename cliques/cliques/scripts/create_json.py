import getopt, sys
import simplejson
from collections import defaultdict
from ipdb import set_trace
import matplotlib.pyplot as plt

# Problems
# Can't handle missing files
# output timestamp
# Can't handle multiple basins

def usage():
    print """-i input prefix e.g. ~/path_graph (where files ~/path_graph_coords.mat etc have been generated
    or you can be specific (will override prefix):
    -e edges_file
    -r energy_file (e.g. stability)
    -c coordinates
    -g graph_file
    -p partitions
    -o output_file
    -b basin_file
    -i input_data (if you want to modify an existing json file)
    """

def parse_args():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "y:ti:b:g:p:r:c:e:h:o:x:", ["help", "output="])
    except getopt.GetoptError, err:
        print str(err)
        usage()
        sys.exit(2)
    
    coords_file = None
    edges_file = None
    energy_file = None
    partitions_file = None
    graph_file = None
    basins_file = None
    timestamp_output = False
    input_file = None
    landscape_type = 'partition'
    output_file = "out"
      
    for o, a in opts:
        if o in ("-t", "--timestamp"):
            import datetime
            timestamp_output = True
            current_time = datetime.datetime.now().isoformat()
            output_file = 'out_{}.json'.format(current_time)
            
    output_file = output_file + ".json"
    
    for o, a in opts:
        if o in ("-x", "--prefix"):
            print "araraad", a
            coords_file = '%s_%s.mat'%(a, 'coords')
            edges_file = '%s_%s.edj'%(a, 'landscape_edgelist')
            energy_file = '%s_%s.mat'%(a, 'energy')
            output_file = '%s_%s.mat'%(a, 'out.json')
            basins_file = '%s_%s.mat'%(a, 'greedy_basins')
            partitions_file = '%s_%s.mat'%(a, 'partitions')
            graph_file = '%s_%s.edj'%(a, 'graph_edgelist')
    
    for o, a in opts:
        if o in ("-c", "--coordinates"):
            coords_file = a
        elif o in ("-e", "--edges"):
            edges_file = a
        elif o in ("-r", "--energy"):
            energy_file = a
        elif o in ("-b", "--basins"):
            basins_file = a
        elif o in ("-o", "--output"):
            output_file = a
        elif o in ("-p", "--partitions"):
            partitions_file = a
        elif o in ("-g", "--graph"):
            graph_file = a
        elif o in ("-i", "--input"):
            input_file = a
        elif o in ("-y", "--type"):
            landscape_type = a
        elif o in ("-h", "--help"):
            usage()
            sys.exit()
   
    return coords_file, edges_file, energy_file, output_file, partitions_file, graph_file, basins_file, input_file, landscape_type

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
    """"""
    f = open(filename, 'r')
    data = []
    for line in f:
        if line:
            l = [cast_type(x) for x in line.split() if x]
            data_point = {'time':l[0], 'values':l[1:]}
            data.append(data_point)
            
    f.close()
    return {'type':'stability', 'data':data}
    
def file_to_basin_process(filename, cast_type):
    """"""
    f = open(filename, 'r')
    data = []
    time_to_values = defaultdict(list)
    for line in f:
        if line:
            l = [x for x in line.strip().split(" ") if x]
            
            nodes = []
            metas = []
            #set_trace()
            try:
                for index, val in enumerate(l[2:]):
                    if index % 2 == 0:
                        nodes.append(int(val))
                    else:
                        metas.append(float(val))
            except:
                print "error"
                #set_trace()            
            time_to_values[float(l[0])].append({'basin':l[1],'nodes':nodes,'metas':metas})
        
    for time in sorted(time_to_values.keys()):
        data.append({'time':float(time), 'values':time_to_values[time]})
            
    f.close()
    return {'type':'basins', 'data':data}

def file_to_basin_sizes(filename, energies, time_limit):
    """"""
    f = open(filename, 'r')
    data = []
     # a dictionary with default value dictionary with default value list
    basin_to_values = defaultdict(lambda: defaultdict(list))
    time_index = 0
    for line in f:
        if line:
            l = [x for x in line.strip().split(" ") if x]
            if float(l[0]) < time_limit:
                #set_trace()
                if float(l[0]) != energies['data'][time_index]['time']:
                    time_index += 1
                    
                basin_size = 0
                basin_prob_size = 0.0
                for index, val in enumerate(l[2:]):
                    if index % 2 == 0:
                        nodes_id = int(val)
                        basin_size += 1
                    else:
                        energy = energies['data'][time_index]['values'][nodes_id]
                        node_to_basin_prob = float(val)
                        basin_prob_size += node_to_basin_prob# * energy
                                              
                basin_to_values[int(l[1])]['time'].append(float(l[0]))
                basin_to_values[int(l[1])]['size'].append(basin_prob_size)
    
    f.close()
    return basin_to_values
    
def main():
    coords_file, edges_file, energy_file, output_file, partitions_file, graph_file, basins_file, input_file, landscape_type = parse_args()
    
    if input_file:
        print "we have input coords"
        f = open(input_file, 'r')
        output = simplejson.load(f)
        f.close()
    else:
        output = defaultdict(list)
    
    if coords_file:
        print "processing coords"
        try:
            output['coords'] = file_to_nested_list(coords_file, float)
        except:
            pass
    if edges_file:
        output['edges'] = file_to_nested_list(edges_file, int)
    if energy_file:
        print "processing energy"
        process = file_to_energy_process(energy_file, float)
        output['processes'].append(process)
    if partitions_file:
        try:
            output['partitions'] = file_to_nested_list(partitions_file, int)
        except:
            print "couldnt open partitions"
    if graph_file:
        import networkx as nx
        graph = nx.read_edgelist(graph_file, nodetype=int, data=(('weight',float),))
        pos = nx.spring_layout(graph,iterations=100)
        output['graph'] = {}
        output['graph']['coords'] = [x.tolist() for x in pos.values()]
        output['graph']['edges'] = file_to_nested_list(graph_file, float)
        plt.figure()
        nx.draw(graph)
    if basins_file:
        process = file_to_basin_process(basins_file, float)
        output['processes'].append(process)

    if landscape_type:
        output['landscapeType'] = landscape_type
    f = open(output_file, 'w')
    simplejson.dump(output, f)

if __name__ == "__main__":
    main()
