// File: main_community.cpp
// -- community detection, sample main file
//-----------------------------------------------------------------------------
// Community detection 
// Based on the article "Fast unfolding of community hierarchies in large networks"
// Copyright (C) 2008 V. Blondel, J.-L. Guillaume, R. Lambiotte, E. Lefebvre
//
// This program must not be distributed without agreement of the above mentionned authors.
//-----------------------------------------------------------------------------
// Author   : E. Lefebvre, adapted by J.-L. Guillaume
// Email    : jean-loup.guillaume@lip6.fr
// Location : Paris, France
// Time	    : February 2008
//-----------------------------------------------------------------------------
// see readme.txt for more details

#include <stdlib.h>
#include <math.h>
#include <string>
#include <iostream> 
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <sys/time.h>


#include "graph_binary.h"
#include "community.h"
extern "C" {
#include <math.h>
#include "mex.h"
#include "matrix.h"
}


using namespace std;

double *data = NULL;
//char *filename_w = NULL;
//char *outfilename = NULL;
int type       = WEIGHTED;
int nb_pass    = 0;
double precision = 0.000001;
int display_level = -1;
int k1 = 16;

void
usage(char *prog_name, char *more) {
  cerr << more;
  cerr << "usage: " << prog_name << " input_file [options]" << endl << endl;
  cerr << "input_file: read the graph to partition from this file." << endl;
  cerr << "-w\t read the graph as a weighted one (weights are set to 1 otherwise)." << endl;
  cerr << "-q epsilon\t a given pass stops when the modularity is increased by less than epsilon." << endl;
  cerr << "-l k\t displays the graph of level k rather than the hierachical structure." << endl;
  cerr << "-h\tshow this usage message." << endl;
  exit(0);
}

/*void
parse_args(int argc, char **argv) {
  if (argc<2)
    usage(argv[0], "Bad arguments number\n");

  for (int i = 1; i < argc; i++) {
    if(argv[i][0] == '-') {
      switch(argv[i][1]) {
      case 'w':
	type = WEIGHTED;
        filename_w = argv[i+1];
	i++;
	break;
      case 'q':
	precision = atof(argv[i+1]);
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
}*/

void
display_time(const char *str) {
  time_t rawtime;
  time ( &rawtime );
  cerr << str << " : " << ctime (&rawtime);
}


extern "C" {
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
 struct timeval tv;
    gettimeofday(&tv,NULL);

  //vector<unsigned long * > pointdegrees;
  //vector<unsigned int * > pointlinks;
  //vector<float * > pointweights;
  srand(((tv.tv_sec*1000)+(tv.tv_usec/1000))*getpid());
  //parse_args(argc, argv);
  //Declarations
    //const mxArray *xData;
    
    
   
    //Copy input pointer x
    //xData = prhs[0];

    //Get matrix x
    data = (double *)mxGetPr(prhs[0]);
    int length_data = mxGetM(prhs[0]);
    
/*
  	// outfilename affectation

	//Declarations
	const mxArray *outfilenameData;
	int outfilenameLength=0;

	//Copy input pointer y
	outfilenameData = prhs[1];

	//Make "TheString" point to the string
	outfilenameLength = mxGetN(outfilenameData)+1;
	outfilename = (char*) mxCalloc(outfilenameLength, sizeof(char)); //mxCalloc is similar to malloc in C
	mxGetString(outfilenameData,outfilename,outfilenameLength);
*/
    
  //time_t time_begin, time_end;
  //time(&time_begin);
  //display_time("start");
    
  vector<vector<int> > output;

  Community *c = new Community(data, length_data, -1, precision);
  int numberinitial=(*c).g.nb_nodes;

  //display_time("file read");

  double mod = (*c).modularity();

  /*cerr << "network : "
       << c.g.nb_nodes << " nodes, " 
       << c.g.nb_links << " links, "
       << c.g.total_weight << " weight." << endl;
 */

  bool improvement = (*c).one_level();

  double new_mod = (*c).modularity();

  //display_time("communities computed");
  //cerr << "modularity increased from " << mod << " to " << new_mod << endl;

  //if (display_level==-1)
    //(*c).display_partition();

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
    c = new Community(g, -1, precision);

	numberofcommunities=(*c).g.nb_nodes;

    /*cerr << "\nnetwork : "
	 << (*c).g.nb_nodes << " nodes, " 
	 << (*c).g.nb_links << " links, "
	 << (*c).g.total_weight << " weight." << endl;*/

    improvement = (*c).one_level();
    new_mod = (*c).modularity();
    
    /*display_time("communities computed");
    cerr << "modularity increased from " << mod << " to " << new_mod << endl;
    */
    //if (display_level==-1)
      //(*c).display_partition();
    output = (*c).display_partition2(output);
    
    g = (*c).partition2graph_binary();
    level++;
    
   // if (level==display_level)
     // g.display();
    

    //display_time("network of communities computed");
  }

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
  //cerr << precision << " " << new_mod << " " << (time_end-time_begin) << endl;

	//Allocate memory and assign output pointer
  
  //free((*c).g.links);
  //free((*c).g.weights);
  //free((*c).g.degrees);
	plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);  //mxReal is our data-type

	//Get a pointer to the data space in our newly allocated memory
	double * out1 = (double*) mxGetPr(plhs[0]);

	out1[0] = new_mod;
    
    
	plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);  //mxReal is our data-type

	//Get a pointer to the data space in our newly allocated memory
	double * out2 = (double*) mxGetPr(plhs[1]);

	out2[0] = numberofcommunities;
    
    int length_output = 0;
    for(int unsigned i=0;i<output.size();i++){
        length_output += output[i].size();
    }


    
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
    delete c;
}
}
