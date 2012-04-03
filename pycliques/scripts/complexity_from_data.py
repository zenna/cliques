from pymatlab.matlab import MatlabSession
from cliques.config import UNIPROJ_HOME
import pickle
import numpy as np

class MyDict(dict):
    def __missing__(self,key):
        self[key] = rv = []
        return rv
results = pickle.load(open('%s/data/karate_range.pickle' % UNIPROJ_HOME,'r'))
x = range(2,18)
partitions_curves = MyDict()
time_curves_x = []
time_curves_y = []
for sample in results:
    for data in sample:
        partitions_curves[data[0]].append(data[1])
        #time_curves_x.append(data[1])
        #time_curves_y.append(float(data[2]))

final_curve = []
for curves in partitions_curves.values():
    final_curve.append(np.mean(curves))
    
session = MatlabSession()
session.run("UNIPROJ_HOME = getenv('UNIPROJ_HOME')")
session.run("save /home/zenna/repos/uniproj/data/karate_2_13_complexity.mat")

session.putvalue('x',np.array(x,dtype='Float64'))
session.putvalue('y',np.array(final_curve,dtype='Float64'))
session.run('semilogy(x,y)')
session.run('hold on')

#session.putvalue('time_x',np.array(time_curves_x,dtype='Float64'))
#session.putvalue('time_y',np.array(time_curves_y,dtype='Float64'))
#session.run('save test.mat')
#todo BELL NUMBERS, 2(n-1), number partitions vs time per partition