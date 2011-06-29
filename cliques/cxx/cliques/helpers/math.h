#pragma once

#include <iostream>
#include <vector>


extern "C" void  dgpadm_ ( int* ideg, int* m, double* t, double* A,
    int* ldh, double* wsp, int* lwsp, int* iwsp,  int* iexp, int* ns, int* iflag);

namespace cliques {

std::vector<double> exp(std::vector<double> matrix, double t, int order) {
    int ideg = 6;
    int m = order;

    int lda=order;
    int ldh=lda;
    int lwsp= 4*ldh*ldh+ideg+1;
    double *wsp = new double[lwsp];
    double *A = &matrix.front();
    int *iwsp = new int[ldh];

    // output
    int iexp, ns, iflag;
    dgpadm_(&ideg, &m, &t, A, &lda, wsp, &lwsp, iwsp, &iexp, &ns, &iflag);

    delete [] wsp;
    delete [] iwsp;

    std::vector<double> output;

    double *start = wsp+iexp-1;


    for (unsigned int i =0;i<matrix.size();++i) {
        output.push_back(*(start+i));
    }
    return output;
}

}
