/*
 * main_varinfo.cpp
 *
 *  Created on: 15 Jan 2010
 *      Author: Antoine Delmotte
 */

#include <stdlib.h>
#include <math.h>
#include <iostream>

extern "C" {
#include "mex.h"
#include "matrix.h"
}

using namespace std;

void quickSort(double arr[], int left, int right) {

	int i = left, j = right;

	double tmp;

	double pivot = arr[(left + right) / 2];

	/* partition */

	while (i <= j) {

		while (arr[i] < pivot)

			i++;

		while (arr[j] > pivot)

			j--;

		if (i <= j) {

			tmp = arr[i];

			arr[i] = arr[j];

			arr[j] = tmp;

			i++;

			j--;

		}

	};

	/* recursion */

	if (left < j)

		quickSort(arr, left, j);

	if (i < right)

		quickSort(arr, i, right);

}

double Max(double *Numbers, int Count)

{

	double Maximum = Numbers[0];

	for (int i = 0; i < Count; i++)

		if (Maximum < Numbers[i])

			Maximum = Numbers[i];

	return Maximum;

}

extern "C" {
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	srand(time(NULL));
	//parse_args(argc, argv);
	//Declarations
	//const mxArray *xData;

	//Copy input pointer x
	//xData = prhs[0];

	//Get matrix x
	double * lnk = (double *) mxGetPr(prhs[0]);
	int nNodes = mxGetM(prhs[0]);
	int nTrials = mxGetN(prhs[0]);
	int *nCommunities = new int[nTrials];

	double *vi_saves;

	if (nlhs > 1) {
		vi_saves = new double[(nTrials) * (nTrials - 1) / 2];
	}

	for (int i = 0; i < nTrials; i++) {
		nCommunities[i] = int(Max(&lnk[i * nNodes], nNodes));
		//printf("ncommunities %d = %d\n", i, nCommunities[i]);
	}

	double infoM = 0;
	double normaM = 0;

	int count = 0;

	for (int x1 = 0; x1 < nTrials; x1++) {
		for (int x2 = 0; x2 < x1; x2++) {

			double * * matrix = new double *[int(nCommunities[x1]) + 1];
			for (int i = 0; i < int(nCommunities[x1]) + 1; i++)
				matrix[i] = new double[int(nCommunities[x2]) + 1];

			double * mh = new double[int(nCommunities[x1]) + 1];
			double * mv = new double[int(nCommunities[x2]) + 1];

			double norma = 0;

			for (int i = 0; i < nCommunities[x1] + 1; i++) {
				mh[i] = 0;
				for (int j = 0; j < nCommunities[x2] + 1; j++) {
					matrix[i][j] = 0;
				}
			}

			for (int j = 0; j < nCommunities[x2] + 1; j++)
				mv[j] = 0;

			for (int i = 0; i < nNodes; i++) {
				matrix[int(lnk[x1 * nNodes + i])][int(lnk[x2 * nNodes + i])]++;
			}

			//printf("matrix = \n");
			for (int i = 0; i < nCommunities[x1] + 1; i++) {
				for (int j = 0; j < nCommunities[x2] + 1; j++) {
					//printf("%f\t", matrix[i][j]);
					mh[i] += matrix[i][j];
					mv[j] += matrix[i][j];
					norma += matrix[i][j];
				}
				//printf("\n");
			}

			//printf("mh = \n");
			//for(int i=0;i<nCommunities[x1]+1;i++)
			//{
			//    printf("%f\t", mh[i]);
			//}

			//printf("mv = \n");
			// for(int i=0;i<nCommunities[x2]+1;i++)
			//{
			//    printf("%f\t", mv[i]);
			//}

			//System.out.println("norma"+norma);

			//for(int i=0;i<nCommunities[oo];i++)
			//	System.out.println(mh[i]);

			double info = 0.0;
			for (int i = 0; i < nCommunities[x1] + 1; i++) {
				for (int j = 0; j < nCommunities[x2] + 1; j++) {
					if (matrix[i][j] > 0)
						info -= 2.0 * matrix[i][j] * log(matrix[i][j] * norma
								/ (mh[i] * mv[j]));
				}
			}

			//printf("info = %f\n", info);

			for (int i = 0; i < nCommunities[x1] + 1; i++)
				free(matrix[i]);
			free(matrix);

			//System.out.println(info);

			double infoh = 0;
			double infov = 0.0;

			// printf("mh_INFO = \n");
			for (int i = 0; i < nCommunities[x1] + 1; i++) {
				infoh += mh[i] * log(mh[i] / norma);
			}

			//printf("%f\n", infoh);

			// printf("mv_INFO = \n");
			for (int j = 0; j < nCommunities[x2] + 1; j++) {
				infov += mv[j] * log(mv[j] / norma);
			}
			// printf("%f\n", infov);

			//printf("norma=%f\n", norma);

			delete[] mh;
			delete[] mv;

			info -= (infoh + infov);
			info /= norma;
			info /= log((double) (nNodes));

			if (info < 0)
				info = 0;

			infoM += info;
			normaM += 1;

			if (nlhs > 1) {
				vi_saves[count] = info;
				count = count + 1;
			}
		}
	}
	infoM /= normaM;

	//printf("infoM = %f\n", infoM);
	delete[] nCommunities;

	//System.out.println(Math.pow(factor,oo-symmetry)+" "+infoM+" "+normaM);

	//System.out.println();


	plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL); //mxReal is our data-type
	//Get a pointer to the data space in our newly allocated memory
	double * out1 = (double*) mxGetPr(plhs[0]);
	out1[0] = infoM;

	if (nlhs > 1) {
		double var_vi = 0;
		for (int i = 0; i < (nTrials * (nTrials - 1)) / 2; i++)
			var_vi += (vi_saves[i] - infoM) * (vi_saves[i] - infoM);
		var_vi = var_vi / (((nTrials) * (nTrials - 1)) / 2 - 1);
		var_vi = sqrt(var_vi);
		plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
		double * out2 = (double*) mxGetPr(plhs[1]);
		out2[0] = var_vi;
	}

	if (nlhs > 5) {
		plhs[5] = mxCreateDoubleMatrix(nTrials, nTrials, mxREAL);
		double * out6 = (double*) mxGetPr(plhs[5]);
		count = 0;
		for (int i = 0; i < nTrials; i++) {
			for (int j = 0; j < nTrials; j++) {
				out6[i * nTrials + j] = 0;
			}
		}
		for (int i = 0; i < nTrials; i++) {
			for (int j = 0; j < nTrials; j++) {
				if (j < i) {
					out6[i * nTrials + j] = vi_saves[count];
					out6[j * nTrials + i] = vi_saves[count];
					count++;
				}
			}
		}
	}

	double median = 0;

	if (nlhs > 2) {
		quickSort(vi_saves, 0, ((nTrials) * (nTrials - 1)) / 2 - 1);
		if (((nTrials) * (nTrials - 1)) / 2 % 2) {
			median = vi_saves[int(ceil(((nTrials) * (nTrials - 1)) / 4.0)) - 1];
		} else {
			median
					= (vi_saves[int(ceil(((nTrials) * (nTrials - 1)) / 4.0))]
							+ vi_saves[int(ceil(((nTrials) * (nTrials - 1))
									/ 4.0)) - 1]) / 2;
		}
		plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL);
		double * out3 = (double*) mxGetPr(plhs[2]);
		out3[0] = median;
	}

	if (nlhs > 3) {
		plhs[3] = mxCreateDoubleMatrix(1, 1, mxREAL);
		double * out4 = (double*) mxGetPr(plhs[3]);
		out4[0] = vi_saves[0];
	}

	if (nlhs > 4) {
		plhs[4] = mxCreateDoubleMatrix(1, 1, mxREAL);
		double * out5 = (double*) mxGetPr(plhs[4]);
		out5[0] = vi_saves[int(((nTrials) * (nTrials - 1)) / 2 - 1)];
	}

	if (nlhs > 6) {
		plhs[6] = mxCreateDoubleMatrix(nTrials * (nTrials - 1) / 2, 1, mxREAL);
		double * out7 = (double*) mxGetPr(plhs[6]);
		count = 0;
		for (int i = 0; i < nTrials * (nTrials - 1) / 2; i++) {
			out7[i] = vi_saves[i];
			count++;
		}

	}

	if (nlhs > 1)
		delete[] vi_saves;

}
}
