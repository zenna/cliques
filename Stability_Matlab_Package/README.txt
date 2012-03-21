###############################################################################
Copyright (C) 2012, <authors> TODO

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.

###############################################################################

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

is included as well as the "Ring-of-rings" graph from the latter publication.
Further example graphs that have been used in these analyses are available on 
request. If you make use of any part of this toolbox, please refer to the 
respective articles.
For detailed instructions on how to compile the code in MATLAB see the file
INSTALL.txt.

If you find a bug or have further comments, please send an email and if 
necessary the input file and the parameters that caused the error.

-----------------------------------------------------------------------------
Authors   : M. Schaub and A. Delmotte
Email     : michael.schaub09@imperial.ac.uk , antoine.delmotte09@imperial.ac.uk 
-----------------------------------------------------------------------------

###############################################################################

-----------------------------------------------------------------------------
Contributions to the code
-----------------------------------------------------------------------------

The C++ code performing the stability optimization is based on the 
implementation of the Louvain method as available from 
http://sites.google.com/site/findcommunities/ 
(Author: Jean-Loup Guillaume)

The code has then been further adapted and extended by R. Lambiotte 
(http://www.lambiotte.be) to allow for the optimization of the stability quality 
function and subsequently been refined by Yun William Yu and Antoine 
Delmotte. 

The MATLAB frontend has been added by Antoine Delmotte. 
Final adjustments and additions, testing, and mainenance is due to 
Antoine Delmotte and Michael Schaub.



