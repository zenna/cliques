import sys
f = open(sys.argv[1])
lines = f.readlines()
for line in lines:
    sl = line.split()
    node1 =  int(sl[0]) -1 
    node2 = int(sl[1]) -1
    weight = int(float(sl[2]) * 10000)
    print node1, node2, weight
    print node2, node1, weight
