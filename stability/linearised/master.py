#!/usr/bin/env python

# multilouvain_final_incr.py
# Author: Yun William Yu <yun.yu09@imperial.ac.uk)
# Date: 2010-03-09
# Purpose: Ensemble maximization of Louvain's clustering stability
# Inputs: 1=graph, 2=total nodes
# Outputs: 
# Copyright (c) 2010 All rights reserved.

# You need Pyro (pyro.sourceforge.net) to run this.
#
# 1) Type "ns" in a shell window to start the Pyro name server.
# 2) Type "python master.py" in a second shell window to start
#    the master process.
# 3) Type "task_manager slave commrunXXX" (where XXX is an
#    integer identifying the process) in a third shell window
#    to start one slave process.

from Scientific.DistributedComputing.MasterSlave import initializeMasterProcess, TaskRaisedException 

#import time

import sys, string, math, subprocess, signal, os
from random import randrange

def signal_handler(signal, frame):
	print 'You pressed Ctrl+C!'
	task.shutdown()
	out_file.close()
	stdout_file.close()
	sys.exit(0)
signal.signal(signal.SIGINT, signal_handler)

pyroid = 'commrun%(pid)03d' %{ 'pid': randrange(1,1000)}
tasks = initializeMasterProcess(pyroid, slave_script="slave.py", use_name_server=True)
job_queue = 24 		# Size of the job queue
numOfSlaves = 4 	# Start off by running this many slave processes (can be added to using task_manager)
tasks.launchSlaveJobs(numOfSlaves)
#for i in range(numOfSlaves):
#	subprocess.Popen(['task_manager','slave',pyroid],stdout=None,stdin=None,stderr=None)

def multilouvain(graph,nodes):
	type(nodes) is int
	print 'Pyro-NS job ID is: ' + pyroid
	print 'multilouvain: \ntime stability Nc starttime counter'
	print  graph, nodes
	out_file = open(graph + '_louvain.out', mode='a')	# If starting a new run, mode='w' for both files
	stdout_file = open(graph + '_louvain.stdout', mode='a')	# Elif resuming a run, mode='a'
	
	for REG in [1,2]:			# REG==1 --> t_m > 1; REG==2 --> t_m < 1
		Ne,Nc = 2,0
		resolution = 300		# Number of timesteps between time t_0 and 10*t_0
		multiplier = math.log(10)/resolution
		u = 0				# Counter (counts number of time steps
		if REG == 1:			# Choose the starting time for t_m > 1
			a = 1 			# Starting time
			Nclim = 1		# Chooses endpoint for clusters found when t_m > 1
		elif REG == 2:			# Choose the starting time for t_m < 1
			a = 1			
			Nclim = nodes		# Chooses endpoint for clusters found when t_m < 1
		else:
			print 'error REG'
			exit(1)
		Rt = Nt = 0
		Nclim_target = 9		# Only end when we've reached Nclim 10 times in a row
		Nclim_count = 0
		#while int(Nc) != int(Nclim):
		while int(Nclim_count) < int(Nclim_target):
			for i in range(job_queue):
				if REG == 1:
					tin_f = a*math.exp(multiplier*u)
					tin = str("%(time)013.7f" % {'time' : tin_f})		# Note that we are storing the time as a string! (for display purposes)
				else:
					tin_f = a*math.exp(-multiplier*u)
					tin = str("%(time)013.7f" % {'time' : tin_f})		# Note that we are storing the time as a string! (for display purposes)
				task_id = tasks.requestTask("communityrun", tin, graph, a, u)
				u = u + 1
			for i in range(job_queue):
				try:
					task_id, tag, result = tasks.retrieveResult("communityrun")
					#print result
					echoline_out = result[0]
					echoline_stdout = result[1]
					Nc = result[2]
					out_file.write(echoline_out+'\n')
					stdout_file.write(echoline_stdout+'\n')
					print echoline_stdout
					if int(Nc) == int(Nclim):
						Nclim_count = Nclim_count + 1	# Add to counter
					else:
						Nclim_count = 0		# Reset
				except TaskRaisedException, e:
					print "Task %s raised %s" % (e.task_id, str(e.exception))
					print e.traceback
			out_file.flush()
			stdout_file.flush()
		#	os.fsync()
	out_file.close()
	stdout_file.close()
	sys.exit(0)

#graph = sys.argv[1]
#nodes = sys.argv[2]
#t0 = time.time()
multilouvain(sys.argv[1], sys.argv[2]) # (4AKE_avg.bin, 3341)
#multilouvain('1MQL2.bin', '22479')
#print time.time() - t0, "seconds"
