/*extern "C" void  dgexpv_ ( int* n, int* m, double* t, double* v, 
    double* w, double* tol, double* anorm,  double* wsp, int* lwsp, int* iwsp, int* liwsp, 
    void (*matvec) ( double*, double*), int* itrace, int* iflag); 
*/

#include <iostream>

extern "C" void  dgpadm_ ( int* ideg, int* m, double* t, double* A, 
    int* ldh, double* wsp, int* lwsp, int* iwsp,  int* iexp, int* ns, int* iflag);

template <typename T>
void print_matrix(T *matrix, int lda) {
    for (int i=0; i<lda; ++i) {
        for (int j=0; j<lda; ++j) {
            std::cout << matrix[j+i*lda] << "\t";
        }
        std::cout << std::endl;
    } 
}

int main() {
    int ideg = 6;
    int m = 10;
    double t = 1.0;
    
    int lda=10;
    int ldh=lda;
    int lwsp= 4*ldh*ldh+ideg+1;
    double *wsp = new double[lwsp];
    int *iwsp = new int[ldh];
    
    // output
    int iexp, ns, iflag;
    double *A = new double[lda*lda];
    
    for (int i =0; i < lda*lda; ++i) {
        if (i%(lda+1) == 0) {
            A[i] = 1.0;
        }
        else {
            A[i] = -1.0;
        }
    }
    
    std::cout << "Input" << std::endl;
    print_matrix(A,lda);
    dgpadm_(&ideg, &m, &t, A, &lda, wsp, &lwsp, iwsp, &iexp, &ns, &iflag);
    
    std::cout << "Input" << std::endl;
    
    print_matrix(wsp+iexp-1,m);
   
    delete [] wsp;
    delete [] A;
    delete [] iwsp;
          
    return 0;
    
}
