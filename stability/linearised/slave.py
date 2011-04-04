# Example for distributed computing using a master-slave setup.
# You need Pyro (pyro.sourceforge.net) to run this example.
#
# 1) Type "ns" in a shell window to start the Pyro name server.
# 2) Type "python master.py" in a second shell window to start
#    the master process.
# 3) Type "task_manager slave demo" in a third shell window
#    to start one slave process.
#
# You can run as many slaves as you want (though for this trivial example,
# the first slave will do all the work before you have time to start a
# second one), and you can run them on any machine on the same local
# network as the one that runs the master process.
#
# See the Pyro manual for other setups, e.g. running slaves on remote
# machines connected to the Internet.
#
# Also see master_slave_demo.py to see how both master and slave can be
# combined within a single script, which is more convenient for short
# scripts.
#

from Scientific.DistributedComputing.MasterSlave import startSlaveProcess
#from Scientific import N
#import subprocess, string
from subprocess import Popen, PIPE
from string import split

# Define (or import) all the task handlers.
def do_communityrun(tin, graph, a, u):
	Re = -10000000
	for i in range(1):	# Loops the community program to find a good partition maximizing stability
		community = Popen(['./community', graph, '-l', '-1', '-t', tin], stdout=PIPE, stderr=PIPE , stdin=None)
		(commT_dat, commMsg) = community.communicate()
		splitted = split(commMsg)
		commE_dat = commT_dat 
		(tout, N, Re, Nc, seconds) = splitted
		I = i
	echoline_out = str(tin)+' '+str(I)+' '+str(Re)+' '+str(Nc)+' '+str(seconds)+' '+str(a)+' '+str(u)
	#echoline_out = '%(tin)s +str(I)+' '+str(Re)+' '+str(Nc)+' '+str(seconds)+' '+str(a)+' '+str(u)
	#out_file.write(echoline_out+'\n)
	hierarchy = Popen(['./hierarchy'],stdin=PIPE,stdout=PIPE)
	hierE_dat = hierarchy.communicate(commE_dat)[0]
	louvain_file = open(graph + '_louvain_' + str(tin) + '.dat', mode='w')
	louvain_file.write(hierE_dat)
	louvain_file.close()
	echoline_stdout = str(tin)+' '+str(Re)+' '+str(Nc)+' '+str(a)+' '+str(u)
	#print echoline_stdout
	#stdout_file.write(echoline+'\n')		
	return ( echoline_out, echoline_stdout, Nc)
		


# Start the slave process after all task handlers have been defined.
startSlaveProcess()
