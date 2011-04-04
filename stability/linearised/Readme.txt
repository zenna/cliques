
These files implement the calculation of the stability and of the 
partitioning with linearised normalised Laplacian over all time scales.

********************************************************************

	INSTALL

********************************************************************

To install the programs, you need a C++ compiler and java compiler  
(should be pre-installed on your mac or linux).

- To install the optimisation code, simply type ./make in your  
terminal (when you're in the directory where the source files are  
located).

- To install the java program, type: javac varinformationdiag.java



********************************************************************

	FILES

********************************************************************

main_convert.cpp : C++ file to be used to convert the graph written
in a txt file into a binary file to be used as an input for
main_community.

The graphs must have the format:
node1 node2 weight_of_bond1
node1 node3 weight_of_bond2
...
when there are links between nodes in the left column and in the right  
column. Important remark: The program does not support self-loops!

./convert -i edgelist.dat -o graph.bin -w

where 

edgelist.dat is the file storing your graph

graph.bin is the binary file of the graph to be used as an input for
main_community.cpp

-w option is to be used if the graph is weighted (don't forget this!)

----------------------------------------------------------------------

main_community.cpp : C++ file implementing the calculation of 
the stability at a particular time and the partitioning. It should be
called in the following way:

./community binary_file -l display level -w -t time

binary_file is the binary file containing the graph
-l option is the display level (put -1 to display all the levels)
-w option is to be put if the graph is weighted
-t option is to specify the markov time

Note also that community takes the options listed in community.cfg,
which is where the number of runs should be specified.

----------------------------------------------------------------------

main_hierarchy.cpp : C++ file organising the communities output of 
main_community.cpp in a table containing a number of lines
equal to the number of nodes and with one column containing the #
of the node and the second the community this node belongs to. It
should be called in the following way:

./hierarchy communities_file
	or
./hierarchy < communities_file

where

communities_file contains the output of main_community.cpp, basically
the nodes and the community they belong to at different levels of the
partitionning.

----------------------------------------------------------------------

varinformationdiag.java : A java program to calculate the variation of 
information between the 50 partitions found for each value of time. 
As an output, you will get the value of this variation of information 
for each value of time. Should be called this way:

	java varinformationdiag nnodes 100

----------------------------------------------------------------------

master.py/slave.py: a script to run the communities calculations
over all time scales. The script will run parallel optimisations 
of t over the entire timescale. Edit the last line of the script with
the name of the graph binary and the total number of nodes. Then make
certain that the Pyro nameserver (pyro-ns) is running before typing
'python master.py'.

Note that Scientific Python (which is distinct from Scipy) must be
installed in order for the script to function.

You can edit the number of subprocesses to spawn in the master.py
script near the top.

NB: master.py does not respond well to Ctrl+C. After force breaking
the program, you'll need to kill the process as well with the 'kill'
command. Then, you'll need to restart pyro-ns if you want to rerun
the program, or else there'll be a namespace clash.


