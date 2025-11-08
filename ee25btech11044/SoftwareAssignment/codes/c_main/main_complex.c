#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include "stb_image.h"
#include "stb_image_write.h"
#include "svd_helpers_complex.h"
int main(int argc,char **argv){
    if(argc<4){
        fprintf(stderr, "Usage: %s input_image output_prefix k1 [k2 ...]\n", argv[0]);
        return 1;
    }
    const char *input_file=argv[1];
    const char *out_prefix=argv[2];
    int width, height, channels;
    unsigned char *imageData = stbi_load(input_file, &width, &height, &channels, 1);
    if(!imageData){
        fprintf(stderr, "Error loading image: %s\n", input_file);
        return 2;
    }
    int rows=height, cols=width;
    double complex *A = (double complex *)malloc(sizeof(double complex) * rows * cols);
    for(int i=0;i<rows * cols;i++)
        A[i] = (double)imageData[i] + 0.0 * I; // real data, complex form
    free(imageData);
    for(int arg=3;arg<argc;arg++){
        int k=atoi(argv[arg]);
        if(k<=0||k>cols) continue;
        double *S = (double *)malloc(sizeof(double) * k);
        double complex *U = (double complex *)malloc(sizeof(double complex) * rows * k);
        double complex *V = (double complex *)malloc(sizeof(double complex) * cols * k);
        double complex *A_k = (double complex *)malloc(sizeof(double complex) * rows * cols);
        int status = computeTruncatedSVDComplex(A, rows, cols, k, S, U, V, 1000, 1e-6);
        if(status != 0){
            fprintf(stderr, "SVD failed for k=%d\n", k);
            free(S);free(U);free(V);free(A_k);
            continue;
        }
        reconstructApproximationComplex(S,U,V,rows,cols,k,A_k);
        double err=computeReconstructionErrorComplex(A,A_k,rows,cols);
        printf("k = %d, Frobenius error = %.3f\n", k, err);
        unsigned char *output = (unsigned char *)malloc(rows * cols);
        for(int i=0;i<rows * cols;i++){
            double val=creal(A_k[i]);
            if(val<0) val=0;
            if(val>255) val=255;
            output[i] = (unsigned char)val;
        }
        char outname[256];
        snprintf(outname, sizeof(outname), "%s_k%d.png", out_prefix, k);
        stbi_write_png(outname,cols,rows,1,output,cols);
        free(output); free(S); free(U); free(V); free(A_k);
    }
    free(A);
    printf("Done\n");
    return 0;
}
