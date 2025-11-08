Project Overview
This project implements Truncated Singular Value Decomposition (SVD) for image compression using C language, based on Power Iteration with Deflation.
The algorithm extracts the top-k
k singular values and vectors of an image matrix and reconstructs a compressed approximation of the original image.
It works for both real and complex matrices and supports grayscale images.

Compilation Instructions (Ubuntu)
1.Go to your main code folder
cd SoftwareAssignment/codes/c_main
2.Compile
gcc -O2 -std=c11 main_complex.c ../c_libs/svd_helpers_complex.c -I../c_libs -lm -o complex_svd

3.Run the program
./a.out ../../figs/globe.jpg ../../figs/globe 5 20 50 100
--The first argument → input image (.jpg or .png)
--The second argument → output prefix for saved images
--The remaining numbers → rank values k for truncated SVD
Example outputs:
globe_k5.png
globe_k20.png
globe_k50.png
globe_k100.png

Output and Analysis
The program prints the Frobenius norm error for each k, showing how reconstruction accuracy improves with higher rank.

Algorithm Summary
The algorithm computes:
1.B=ATA
2.Iteratively finds the top-k eigenvectors of B using Power Iteration.
3.Each eigenvalue’s square root gives a singular value of A
4.Left singular vectors are computed as U=Av∥Av∥
5.Each found component is deflated from B to get the next one.
6.Reconstructs Ak=UkΣkVkT
This process is efficient, accurate for large images, and suitable for both real and complex matrices.

Notes
Developed in C for performance and memory control.
Uses stb_image libraries for portable image I/O.
Works on any grayscale image.

Reference
Gilbert Strang, "The Singular Value Decomposition and the Geometry of Matrices", MIT OpenCourseWare.
