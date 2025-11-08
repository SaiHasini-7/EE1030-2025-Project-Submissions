#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "svd_helpers_complex.h"
static void computeHermitianProduct(double complex *A,double complex *B,int m,int n){
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            double complex sum= 0.0 + 0.0 * I;
            for(int k=0;k<m;k++)
                sum += conj(A[k * n + i]) * A[k * n + j];
            B[i * n + j]=sum;
        }
    }
}
static void normalizeVector(double complex *v,int n){
    double norm=0.0;
    for(int i=0;i<n;i++) norm += pow(cabs(v[i]), 2);
    norm=sqrt(norm);
    if(norm==0.0) return;
    for(int i=0;i<n;i++) v[i] /= norm;
}
int computeTruncatedSVDComplex(double complex *A,int m,int n,int k,double *S,double complex *U, double complex *V,int max_iter, double tol){
    double complex *B = (double complex *)malloc(sizeof(double complex) * n * n);
    computeHermitianProduct(A,B,m,n);
    for (int comp=0;comp<k;comp++){
        double complex *v = (double complex *)malloc(sizeof(double complex) * n);
        for(int i=0;i<n;i++)
            v[i] = (double)rand() / RAND_MAX + ((double)rand() / RAND_MAX) * I;
        normalizeVector(v, n);
        for(int iter=0;iter<max_iter;iter++){
            double complex *w = (double complex *)malloc(sizeof(double complex) * n);
            for(int i=0;i<n;i++){
                double complex sum= 0.0 + 0.0 * I;
                for(int j=0;j<n;j++)
                    sum += B[i * n + j] * v[j];
                w[i]=sum;
            }
            normalizeVector(w, n);
            double diff = 0.0;
            for(int i=0;i<n;i++) diff += pow(cabs(w[i] - v[i]), 2);
            diff=sqrt(diff);
            free(v);
            v=w;
            if(diff<tol) break;
        }
        double complex lambda = 0.0 + 0.0 * I;
        for(int i=0;i<n;i++) {
            double complex temp = 0.0 + 0.0 * I;
            for(int j=0;j<n;j++)
                temp += B[i * n + j] * v[j];
            lambda += conj(v[i]) * temp;
        }
        double sigma = sqrt(creal(lambda));
        S[comp]=sigma;

        double complex *u = (double complex *)malloc(sizeof(double complex) * m);
        for(int i=0; i<m;i++){
            double complex sum = 0.0 + 0.0 * I;
            for(int j=0;j<n;j++)
                sum += A[i * n + j] * v[j];
            u[i]=sum;
        }
        normalizeVector(u, m);
        for(int i=0;i<m;i++) U[i * k + comp] = u[i];
        for(int i=0;i<n;i++) V[i * k + comp] = v[i];
        free(u);
        for(int i=0;i<n;i++){
            for(int j=0;j<n;j++){
                B[i * n + j] -= sigma * sigma * v[i] * conj(v[j]);
            }
        }
        free(v);
    }
    free(B);
    return 0;
}
void reconstructApproximationComplex(double *S,double complex *U,double complex *V,int m,int n,int k,double complex *A_reconstructed){
    for(int i=0;i<m;i++){
        for(int j=0;j<n;j++){
            double complex sum = 0.0 + 0.0 * I;
            for(int c=0;c<k;c++)
                sum += S[c] * U[i * k + c] * conj(V[j * k + c]);
            A_reconstructed[i * n + j]=sum;
        }
    }
}
double computeReconstructionErrorComplex(double complex *A,double complex *A_k,int m,int n){
    double err=0.0;
    for(int i=0;i<m * n;i++){
        double diff = cabs(A[i] - A_k[i]);
        err += diff * diff;
    }
    return sqrt(err);
}
