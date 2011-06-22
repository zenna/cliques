import getopt, sys
import simplejson
from collections import defaultdict

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
        opts, args = getopt.getopt(sys.argv[1:], "ti:b:g:p:r:c:e:h:o:x:", ["help", "output="])
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
    input_data = None
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
            coords_file = '{}_{}.mat'.format(a, 'coords')
            edges_file = '{}_{}.mat'.format(a, 'edges')
            energy_file = '{}_{}.mat'.format(a, 'energy')
            output_file = '{}_{}.mat'.format(a, 'out.json')
            basins_file = '{}_{}.mat'.format(a, 'basins')
            partitions_file = '{}_{}.mat'.format(a, 'partitions')
            graph_file = '{}_{}.edj'.format(a, 'graph')
    
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
        elif o in ("-h", "--help"):
            usage()
            sys.exit()
        else:
            assert False, "unhandled option"
   
    return coords_file, edges_file, energy_file, output_file, partitions_file, graph_file, basins_file, input_data
            
def file_to_nested_list(filename, cast_type):
    """Convert file to matrix"""
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
            for index, val in enumerate(l[2:]):
                if index % 2 == 0:
                    nodes.append(int(val))
                else:
                    metas.append(cast_type(val))
                        
            time_to_values[l[0]].append({'basin':l[1],'nodes':nodes,'metas':metas})
        
    for time in sorted(time_to_values.keys()):
        data.append({'time':float(time), 'values':time_to_values[time]})
            
    f.close()
    return {'type':'basins', 'data':data}
    
def main():
    coords_file, edges_file, energy_file, output_file, partitions_file, graph_file, basins_file, input_file = parse_args()
    
    if input_file:
        f = open(input_file, 'r')
        output = json.load(f)
        f.close()
    else:
        output = defaultdict(list)
    
    if coords_file:
        output['coords'] = file_to_nested_list(coords_file, float)
    if edges_file:
        output['edges'] = file_to_nested_list(edges_file, int)
    if energy_file:
        process = file_to_energy_process(energy_file, float)
        output['processes'].append(process)
    if partitions_file:
        output['partitions'] = file_to_nested_list(partitions_file, int)
    if graph_file:
        import networkx as nx
        graph = nx.read_edgelist(graph_file, nodetype=int, data=(('weight',float),))
        pos = nx.spring_layout(graph,iterations=100)
        output['graph'] = {}
        output['graph']['coords'] = [x.tolist() for x in pos.values()]
        output['graph']['edges'] = file_to_nested_list(graph_file, float)
    if basins_file:
        process = file_to_basin_process(basins_file, float)
        output['processes'].append(process)
        
    f = open(output_file, 'w')
    simplejson.dump(output, f)

if __name__ == "__main__":
    main()
