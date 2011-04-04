
#include <stdio.h>
#include <math.h>

extern "C" {
#include <math.h>
#include "mex.h"
#include "matrix.h"
}

//using namespace std;


double maximum(double *data, int lng_data){
    double max=0;
    for(int i=1;i<lng_data;i++){
        if(data[i]>max)
            max=data[i];
    }
    return max;
}

extern "C" {
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
    double * data = (double *)mxGetPr(prhs[0]);
    int nNodes = mxGetM(prhs[0]);
    
    
    int count=0;
    //int nNodes= int(sqrt(length_data));

    double limit=maximum(data,nNodes*nNodes)/1000000;
    
    for(int i=0;i<nNodes;i++){
        for(int j=0;j<nNodes;j++){
            if(data[i*nNodes+j]>limit)
                count++;
        }
    }
    
    plhs[0] = mxCreateDoubleMatrix(count,3,mxREAL);
    double * edges = (double*) mxGetPr(plhs[0]);
    
    int count2=0;
    for(int i=0;i<nNodes;i++){
        for(int j=0;j<nNodes;j++){
             if(data[i*nNodes+j]>limit){
                edges[count2]=i;
                edges[count+count2]=j;
                edges[2*count+count2]=data[i*nNodes+j];
                count2++;
             }
        }
    }
    
}
}

