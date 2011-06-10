import getopt, sys
import simplejson

def usage():
    print """-i input prefix e.g. ~/path_graph (where files ~/path_graph_coords.mat etc have been generated
    or you can be specific (will override prefix):
    -e edges_file
    -r energy_file (e.g. stability)
    -c coordinates
    -o output_file
    """

def parse_args():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "p:r:c:e:h:o:", ["help", "output="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err)
        usage()
        sys.exit(2)
    output_file = 'out.json'
    coords_file = None
    edges_file = None
    energy_file = None
    partitions_file = None
    
    for o, a in opts:
        if o in ("-i", "--input"):
            coords_file = '{}_{}.mat'.format(a, 'coords')
            edges_file = '{}_{}.mat'.format(a, 'edges')
            energy_file = '{}_{}.mat'.format(a, 'energy')
            output_file = '{}_{}.mat'.format(a, 'out.json')
    
    for o, a in opts:
        if o in ("-c", "--coordinates"):
            coords_file = a
        elif o in ("-e", "--edges"):
            edges_file = a
        elif o in ("-r", "--energy"):
            energy_file = a
        elif o in ("-o", "--output"):
            output_file = a
        elif o in ("-p", "--partitions"):
            partitions_file = a
        elif o in ("-h", "--help"):
            usage()
            sys.exit()
        else:
            assert False, "unhandled option"
   
    return coords_file, edges_file, energy_file, output_file, partitions_file
            
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

def main():
    output = {}
    coords_file, edges_file, energy_file, output_file, partitions_file = parse_args()
    if coords_file:
        output['coords'] = file_to_nested_list(coords_file, float)
    if edges_file:
        output['edges'] = file_to_nested_list(edges_file, int)
    if energy_file:
        output['energies'] = file_to_nested_list(energy_file, float)
    if partitions_file:
        output['partitions'] = file_to_nested_list(partitions_file, int)
    
    f = open(output_file, 'w')
    simplejson.dump(output, f)

if __name__ == "__main__":
    main()
