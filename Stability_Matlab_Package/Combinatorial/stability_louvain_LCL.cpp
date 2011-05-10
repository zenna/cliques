// File: main_community.cpp
// -- community detection, sample main file
//-----------------------------------------------------------------------------
// Community detection 
// Based on the article "Fast unfolding of community hierarchies in large networks"
// Copyright (C) 2008 V. Blondel, J.-L. Guillaume, R. Lambiotte, E. Lefebvre
// and on the article "Dynamics and mModular Structure in Networks"
// Copyright (C) 2008 R. Lambiotte, J.-C. Delvenne and M. Barahona
//
// This program must not be distributed without agreement of the above mentionned authors.
//-----------------------------------------------------------------------------
// Author   : E. Lefebvre, adapted by J.-L. Guillaume and R. Lambiotte
// Email    : r.lambiotte@imperial.ac.uk
// Location : London, UK
// Time	    : November 2008
//-----------------------------------------------------------------------------
// see readme.txt for more details

#if defined(macintosh) || defined(MACOS) || defined(_MACOS)  || defined(__APPLE__)  || defined(_MAC)   || defined(MAC)  || defined(mac)  || defined(MACINTOSH)
#define __MAC__
#endif


#if defined(_WIN32)  || defined(WIN32)  || defined(__WIN32__) || defined(_WIN64)  || defined(WIN64) || defined(__WIN64__) || defined(_WINDOWS) || defined(__WINDOWS__)
#define __WIN__
#endif

#if defined(linux)  || defined(__linux__)  || defined(_linux) || defined(LINUX)  || defined(_LINUX) || defined(_UNIX) || defined(__UNIX__) || defined(__gnu_linux__) || defined(__unix__) || defined(UNIX) || defined(unix) || defined(sparc)
#define __lin__
#endif

#include <stdlib.h>
#include <math.h>
#include <string>
#include <iostream> 
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <string.h>
#include <stdio.h>


#ifdef __WIN__
#include <time.h>
#include <sys/timeb.h> 
#include <process.h>
#include <windows.h>
#endif


#include "graph_binary.h"
#include "community.h"
extern "C" {
#include <math.h>
#include "mex.h"
#include "matrix.h"

}

#ifdef __WIN__
int gettimeofday (struct timeval *tp, void *tz) 
{ 
struct _timeb timebuffer; 
_ftime (&timebuffer); 
tp->tv_sec = timebuffer.time; 
tp->tv_usec = timebuffer.millitm * 1000; 
return 0; 
} 
#endif

#ifdef __lin__
#include <sys/time.h>
#endif

#ifdef __MAC__
#include <sys/time.h>
#endif




using namespace std;

double *data = NULL;
//char *filename_w = NULL;
//char *outfilename = NULL;
int type       = WEIGHTED;
int nb_pass    = 0;
double precision = 0.000001;
int display_level = -1;
//int k1 = 16;
int length_data=-1;
double timet = 1.0;
bool hierarchy = false;

/*void
usage(char *prog_name, char *more) {
  cerr << more;
  cerr << "usage: " << prog_name << " input_file [options]" << endl << endl;
  cerr << "input_file: read the graph to partition from this file." << endl;
  cerr << "-w\t read the graph as a weighted one (weights are set to 1 otherwise)." << endl;
  cerr << "-q epsilon\t a given pass stops when the modularity is increased by less than epsilon." << endl;
  cerr << "-t Optimisation of a quantity generalzing modularity in order to uncover smaller communities. Modularity: t=1. t must be smaller than 1." << endl;
  cerr << "-l k\t displays the graph of level k rather than the hierachical structure." << endl;
  cerr << "-h\tshow this usage message." << endl;
  exit(0);
}

void
parse_args(int argc, char **argv) {
  if (argc<2)
    usage(argv[0], "Bad arguments number\n");

  for (int i = 1; i < argc; i++) {
    if(argv[i][0] == '-') {
      switch(argv[i][1]) {
      case 'w':
	type = WEIGHTED;
	break;
      case 'q':
	precision = atof(argv[i+1]);
	i++;
	break;
	  case 't':
	timet = atof(argv[i+1]);
	i++;
	break;
      case 'l':
	display_level = atoi(argv[i+1]);
	i++;
	break;
      case 'k':
	k1 = atoi(argv[i+1]);
	i++;
	break;
      default:
	usage(argv[0], "Unknown option\n");
      }
    } else {
      if (filename==NULL)
        filename = argv[i];
      else
        usage(argv[0], "More than one filename\n");
    }
  }
}

void
display_time(const char *str) {
  time_t rawtime;
  time ( &rawtime );
  cerr << str << " : " << ctime (&rawtime);
}
*/

bool parse_arg(int nrhs, const mxArray *prhs[]){
    
	if(nrhs>0){
        if(mxGetN(prhs[0])!=3 && mxGetN(prhs[0])!=2){
            printf("N=%d",mxGetN(prhs[0]));
            return false;
        }
		data = (double *)mxGetPr(prhs[0]);
  		length_data = mxGetM(prhs[0]);
	}
    if(nrhs>1){
        timet=((double) mxGetScalar(prhs[1]));
	}
    if(nrhs>2){
        if(precision>1)
            return false;
		precision=((double) mxGetScalar(prhs[2]));
		//precision=tempo;
     }
    if(nrhs>3){
        double p = (double) mxGetScalar(prhs[3]);
		if (p==119){
            if(mxGetN(prhs[0])!=3)
                return false;
			type = WEIGHTED;
        }else if(p==117){
            type = UNWEIGHTED;
        }else{
            return false;
        }
    }
    if(nrhs>4){ 
        double p = (double) mxGetScalar(prhs[4]);
		if (p==104){
			hierarchy=true;
        }else if(p==110){
            hierarchy=false;
        }else{
            return false;
        }
	}
    if(nrhs>5 || nrhs<1){
        return false;
	}
    return true;
}

extern "C" {
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
  
struct timeval tv;
gettimeofday(&tv,NULL);


  // Initialize the random generator
  srand(((tv.tv_sec * 1000) + (tv.tv_usec / 1000))*getpid());

  //parse arguments
  if(parse_arg(nrhs, prhs) && nlhs <4){
    

    
  // Vector containing the hierarchy of the partitioning    
  vector<vector<int> > output;

  // Creation of the Community object 
  Community *c = new Community(data, length_data, -1, precision, timet, type);
  

  // Initial number of nodes
  int numberinitial=(*c).g.nb_nodes;

  // Calculation of the initial modularity
  double mod = (*c).modularity();

  // First partitioning
  bool improvement = (*c).one_level();
  
  // Calculation of its modularity
  double new_mod = (*c).modularity();

  //display_time("communities computed");
  //cerr << "modularity increased from " << mod << " to " << new_mod << endl;

  //if (display_level==-1)
    //(*c).display();

  //output.push_back((*c).n2c); 
  output = (*c).display_partition2(output);

  Graph g = (*c).partition2graph_binary();

 
//if (display_level==0)
  //  g.display();

  //display_time("network of communities computed");
  int numberofcommunities=(*c).g.nb_nodes;
  int level=0;
  while(improvement) {
    //pointdegrees.pushback((*c).g.degrees);
    //pointlinks.pushback((*c).g.links);
    //pointweights.pushback((*c).g.weights);
    mod=new_mod;
    
    //free((*c).g.links);
    //free((*c).g.degrees);
    //free((*c).g.weights);

    delete c;
    c = new Community(g, -1, precision, timet, numberinitial);

    
    numberofcommunities=(*c).g.nb_nodes;

    /*cerr << "\nnetwork : "
	 << (*c).g.nb_nodes << " nodes, " 
	 << (*c).g.nb_links << " links, "
	 << (*c).g.total_weight << " weight." << endl;*/
    
    improvement = (*c).one_level();
    
    
    if ((*c).nb_pass!=-10)
        new_mod = (*c).modularity();
    else
        new_mod =0;
    
    
    //new_mod = (*c).modularity();
    
    /*display_time("communities computed");
    cerr << "modularity increased from " << mod << " to " << new_mod << endl;
    */
    //if (display_level==-1)
      //(*c).display();
      
    output = (*c).display_partition2(output);
    
   
    g = (*c).partition2graph_binary();
    level++;

   // if (level==display_level)
     // g.display();
    

    //display_time("network of communities computed");
  }

  numberofcommunities=(*c).g.nb_nodes;
  free(g.links);
  free(g.weights);
  free(g.degrees);
  free(g.nb_nodes_per_comm);
    /*pointdegrees.pushback((*c).g.degrees);
    pointlinks.pushback((*c).g.links);
    pointweights.pushback((*c).g.weights);

	for(int i=0;i<pointdegrees.size();i++){
	if(pointdegrees[i]!=NULL)
		free(pointdegrees[i]);
	}
	for(int i=0;i<pointlinks.size();i++){
	if(pointlinks[i]!=NULL)
		free(pointlinks[i]);
	}
	for(int i=0;i<pointweights.size();i++){
	if(pointweights[i]!=NULL)
		free(pointweights[i]);
	}*/
  //time(&time_end);

  //cerr << 0 << " " << numberinitial << " " << new_mod << " "<< numberofcommunities << " " << (time_end-time_begin) << endl;
//	cout << "" << endl;
  //cerr << precision << " " << new_mod << " " << (time_end-time_begin) << end/home/apd09/Documentsl;

	//Allocate memory and assign output pointer
  
  //free((*c).g.links);
  //free((*c).g.weights);
  //free((*c).g.degrees);
	plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);  //mxReal is our data-type

	//Get a pointer to the data space in our newly allocated memory
	double * out1 = (double*) mxGetPr(plhs[0]);
    
	out1[0] = new_mod;
    
    if(nlhs>1){
    
        plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);  //mxReal is our data-type

        //Get a pointer to the data space in our newly allocated memory
        double * out2 = (double*) mxGetPr(plhs[1]);

        out2[0] = numberofcommunities;


    }


    if(hierarchy && nlhs>2){    
        int length_output = 0;
        for(int unsigned i=0;i<output.size();i++)
            length_output += output[i].size();
        plhs[2] = mxCreateDoubleMatrix(length_output,2,mxREAL);
        double * output_tab = (double*) mxGetPr(plhs[2]);
        double counter_temp = 0;
        int counter_temp2 = 0;
        for(int i=0;i<output.size();i++){
            counter_temp = 0;
            for(int unsigned j=0;j<output[i].size();j++){

                output_tab[counter_temp2] = (double) counter_temp;
                output_tab[counter_temp2+length_output] = (double) output[i][j];
                counter_temp2 ++;
                counter_temp ++;
            }
        }
    }else if(nlhs>2){
        
        vector<int> n2c(numberinitial);
      
        for (unsigned int i=0 ; i<numberinitial ; i++)
          n2c[i]=i;

        for (int l=0 ; l<output.size() ; l++){
              for (unsigned int node=0 ; node<numberinitial ; node++){
                n2c[node] = output[l][n2c[node]];
            }
        }

        plhs[2] = mxCreateDoubleMatrix(numberinitial,1,mxREAL);
        double * output_tab = (double*) mxGetPr(plhs[2]);
        for (unsigned int node=0 ; node<numberinitial ; node++){
            output_tab[node] = (double) n2c[node];

        }    
    }
    
    delete c;
    
  }else{
      printf("\n\nSYNTAX:");
      
      printf("\n\n\t[stability]=stability_louvain_LCL(graph)");
      printf("\n\t[stability]=stability_louvain_LCL(graph, markov time)");
      printf("\n\t[stability]=stability_louvain_LCL(graph,  markov time, precision)");
      printf("\n\t[stability]=stability_louvain_LCL(graph,  markov time, precision, weighted)");
      printf("\n\t[stability]=stability_louvain_LCL(graph,  markov time, precision, weighted, hierarchy)");
      printf("\n\t[stability, number of communities]=stability_louvain_LCL(graph)");
      printf("\n\t[stability, number of communities]=stability_louvain_LCL(graph, markov time)");
      printf("\n\t[stability, number of communities]=stability_louvain_LCL(graph,  markov time, precision)");
      printf("\n\t[stability, number of communities]=stability_louvain_LCL(graph,  markov time, precision, weighted)");
      printf("\n\t[stability, number of communities]=stability_louvain_LCL(graph,  markov time, precision, weighted, hierarchy)");
      printf("\n\t[stability, number of communities, communities]=stability_louvain_LCL(graph)");
      printf("\n\t[stability, number of communities, communities]=stability_louvain_LCL(graph, markov time)");
      printf("\n\t[stability, number of communities, communities]=stability_louvain_LCL(graph,  markov time, precision)");
      printf("\n\t[stability, number of communities, communities]=stability_louvain_LCL(graph,  markov time, precision, weighted)");
      printf("\n\t[stability, number of communities, communities]=stability_louvain_LCL(graph,  markov time, precision, weighted, hierarchy)");
      
      printf("\n\nDESCRIPTION:");
      
      printf("\n\n  Inputs:");
	  printf("\n\n\t - graph\t list of edges in the form [node, node, weight;node, node, weight;...] \n\t\t\t if the graph is weighted or [node, node;node, node;...] if the graph is unweighted.");
      printf("\n\t - markov time\t time allowed to the random walkers to explore the graph \n\t\t\t (default: 1.0, corresponding to the calculation of the modularity)");
      printf("\n\t - precision\t required precision for the stability (default: 0.000001)");
      printf("\n\t - weighted\t put \'w\' if the graph is weighted or \'u\' if the graph is unweighted (default: unweighted)");
      printf("\n\t - hierarchy\t put \'h\' if the ''communities'' output should be the hierarchy of \n\t\t\t the partitions or \'n\' if it should only be the final partition (default: non-hiearchy)");
      
      printf("\n\n  Outputs:");
      printf("\n\n\t - stability\t\t\t value of the stability of the final partition.");
      printf("\n\t - number of communities\t number of communities of the final partition.");
      printf("\n\t - communities\t\t\t affectation of the nodes to the communities either in \n\t\t\t\t\t the final partition (if \'hierarchy\'==\'u\') or in the hierarchy of partitions (if \'hierarchy\'==\'h\').");
      printf("\n\n\n");
      
      plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
      
      if(nlhs>1)
          plhs[1] = mxCreateDoubleMatrix(0,0,mxREAL);
      
      if(nlhs>2)
        plhs[2] = mxCreateDoubleMatrix(0,0,mxREAL);
  }
  
  
  
}
}
