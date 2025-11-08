#ifndef SVD_HELPERS_COMPLEX_H
#define SVD_HELPERS_COMPLEX_H
#include <complex.h>
int computeTruncatedSVDComplex(double complex *A,int m,int n,int k,double *S,double complex *U,double complex *V,int max_iter, double tol);
void reconstructApproximationComplex(double *S,double complex *U,double complex *V,int m,int n,int k,double complex *A_reconstructed);
double computeReconstructionErrorComplex(double complex *A, double complex *A_k,int m, int n);
#endif
