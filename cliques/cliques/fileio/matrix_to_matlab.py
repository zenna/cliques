import sys
f = open(sys.argv[1])
lines = f.readlines()
for i, line in enumerate(lines):
    weights = line.split()
    
    for j, weight in enumerate(weights):
        if i != j:
            print i, j, int(float(weight)*10000)