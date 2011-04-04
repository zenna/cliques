

These files implement the calculation of the stability and of the 
partitioning with full normalised Laplacian over all time scales.

********************************************************************

	INSTALL

********************************************************************

Execute the following lines in Matlab:

First, verify the compiler you are using by typing 
	
	mex -setup

and choose a compiler (I don't think it works with the LCC compiler).

Then, compile the mex files by typing:

	mex main_community_mex.cpp community.cpp graph_binary.cpp
	mex main_varinfo.cpp
	mex main_hierarchy.cpp

********************************************************************

	FILES

********************************************************************

main_community_mex.cpp : Mex file implementing the calculation of 
the modularity (stability at time = 1) and the partitioning. Should 
be called from Matlab in the following way:

	[stability, nb_comm, communities] = main_community_mex(graph);

where 

stability is the modularity of the partitioning,

nb_comm is the number of communities in the partitioning,

communities is a matrix (# nodes + # communities 1st level + 
# communities 2nd level + ..., 2), the first column containing the 
number of the node and the second, the cluster the node belongs to.

graph is a matrix (# bonds*2, 3) containing the list of the bonds 
in the graph. The first two columns contain the number of the nodes 
involved in the considered bonds and the third column contains the 
weight of the bond.

-------------------------------------------------------------------

main_varinfo.cpp : C++ file implementing the calculation of the 
variation of information between the different partitions found. 
It should be called from Matlab in the following way:

	var_info = main_varinfo(lnk);

where 

var_info is the variation of variation between the different 
partitions obtained

lnk is a matrix (# nodes, # trials) containing in each column the
cluster to which each node belongs to (ith line contains the # of the
cluster for node i).

-------------------------------------------------------------------

main_hierarchy.cpp : Mex file organising the communities output of 
main_communities_mex.cpp in a matrix containing a number of lines
equal to the number of nodes and with one column containing the #
of the node and the second the community this node belongs to. It
should be called in the following way:

communities_hierarchy = main_hierarchy(communities);

where 

communities is the output of main_community_mex.cpp

communities_hierarchy is the matrix containing the cluster each node
belongs to.

-------------------------------------------------------------------

script_comm_all_times.m : Matlab script implementing the calculation 
of the communities over all time scales.

-------------------------------------------------------------------

expmat.m : Matlab function implementing the calculation of the 
exponential of the Normalised Laplacian at a particular time. 
This generates a new graph of which modularity is the stability of 
the original graph at this time.

	graph = expmat(data,time);

where

data is the name of the mat file containing the graph (the variable
should be called Graph).

time is the time at which the stability is to be calculated

graph is the graph generated and of which modularity is the stability
time "time".

-------------------------------------------------------------------

expmatprime.m : Matlab function implementing the calculation of the 
exponential of the Standard Laplacian at a particular time. 
This generates a new graph of which modularity is the stability of 
the original graph at this time.

	graph = expmatprime(data,time);

where

data is the name of the mat file containing the graph (the variable
should be called Graph).

time is the time at which the stability is to be calculated

graph is the graph generated and of which modularity is the stability
time "time".

--------------------------------------------------------------------








