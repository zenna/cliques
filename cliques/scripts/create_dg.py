

maxima = [22857,22860,22863,23042,23072,23107,31227,33656,33800,34053,34213,34248,67290,69785,69788,69791,69970,70000,70035,111838,123953,123988,125241,125389,127643,127646,127649,127828,127858,127893,150933,156216,159420,159833,159868,163798]

l_maxima = {}
maxima_values = {}

import numpy as np

i = 0
for maximum in maxima:
    l_maxima[maximum] = i
    file = open('/home/zenna/repos/uniproj/data/stabilities/renaud_n14.stabilities')
    for x in range(maximum+2):
        lines = file.readlines()
        
    A = line.split(',')
    maxima_values[i] = A[124]
    file.close()
    print maxima_values[i]

    i = i + 1

#file = open('/home/zenna/repos/uniproj/data/bottlenecks/temp.bsp')
#lines = file.readlines()


