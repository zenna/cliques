-----------------------------------------------------------------------------
Community Detection using the stability of a graph partition.
-----------------------------------------------------------------------------

The code implements the stability method as discussed in the articles

"Stability of graph communities across time scales" 
Delvenne, J.-C.; Yaliraki, S. N. & Barahona, M.
Proceedings of the National Academy of Sciences, 2010, 107, 12755-12760 

and

"Laplacian Dynamics and Multiscale Modular Structure in Networks"
Lambiotte, R.; Delvenne, J.-C. & Barahona, M.
http://arxiv.org/abs/0812.1770, 2009



For optimizing the stability quality function the Louvain algorithm as described 
in the publication

"Fast unfolding of communities in large networks",
Vincent D Blondel, Jean-Loup Guillaume, Renaud Lambiotte, Etienne Lefebvre,
Journal of Statistical Mechanics: Theory and Experiment 2008 (10), P10008

is used.
A simple example graph is included in the folder /demo/ that demonstrates the 
main functionality of the program.
As a second example a graph representation of Adenylate Kinase (AdK) as it has 
been used in

"Protein multi-scale organization through graph partitioning and robustness 
analysis: application to the myosin–myosin light chain interaction"
Delmotte, A.; Tate, E. W.; Yaliraki, S. N. & Barahona, M. 
Physical Biology, 2011, 8, 055010

and

"Markov dynamics as a zooming lens for multiscale community detection: 
non clique-like communities and the field-of-view limit"
Schaub, M. T.; Delvenne, J.-C.; Yaliraki, S. N. & Barahona, M. 
http://arxiv.org/abs/1109.5593, 2011

Further example graphs that have been used in these analyses are available on 
request. If you make use of any part of this toolbox, please refer to the 
respective articles.
For detailed instructions on how to compile the code in MATLAB see the file
INSTALL.txt.

-----------------------------------------------------------------------------
Authors   : 
Email     : 
-----------------------------------------------------------------------------

###############################################################################

Disclaimer:
This program or any part of it must not be distributed without prior agreement 
of the above mentionned authors.
This is the first public version of this program, although we have made an 
effort to make it compatible with all possible OS, we give no garantuees that it
will work with every possible system setup.

If you find a bug or have further comments, please send an email to TODO
and if necessary the input file and the parameters that caused the error.

###############################################################################
