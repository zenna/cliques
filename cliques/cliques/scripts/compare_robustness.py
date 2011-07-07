import getopt, sys
import simplejson
from collections import defaultdict
import matplotlib.pyplot as plt

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
            edges_file = '{}_{}.edj'.format(a, 'landscape_edgelist')
            energy_file = '{}_{}.mat'.format(a, 'energy')
            output_file = '{}_{}.mat'.format(a, 'out.json')
            basins_file = '{}_{}.mat'.format(a, 'greedy_basins')
            partitions_file = '{}_{}.mat'.format(a, 'partitions')
            graph_file = '{}_{}.edj'.format(a, 'graph_edgelist')
    
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
   
    return coords_file, edges_file, energy_file, output_file, partitions_file, graph_file, basins_file, input_data
            
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
    
def main():
    coords_file, edges_file, energy_file, output_file, partitions_file, graph_file, basins_file, input_file = parse_args()
    
    # Load 
    if basins_file:
        try:
            process = file_to_basin_process(basins_file, float)
            'time':float(time v({'basin':l[1],'nodes':nodes,'metas':metas
            output['processes'].append(process)
        except:
            print "no basin file"


    data = {}
    plt.semilogx([1,10,100,1000], [1,2,3, 4], 'go-', label='variation of information', linewidth=2)
    plt.semilogx([1,10,100], [1,2,5], 'ro-', label='robustness', linewidth=2)
    plt.legend()
    plt.show()
    
if __name__ == "__main__":
    main()
