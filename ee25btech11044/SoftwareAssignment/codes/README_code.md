Description
This code implements Truncated Singular Value Decomposition (SVD) using the
Power Iteration with Deflation method in C

Compilation

Run the following commands in Ubuntu/Linux:
cd SoftwareAssignment/codes/c_main
gcc -O2 -std=c11 main_complex.c ../c_libs/svd_helpers_complex.c -I../c_libs -lm -o complex_svd
-lm -o svd_app
./svd_app <input_image> <output_prefix> <k1> [k2...]

Notes
Works with .png, .jpg, and .jpeg images.
All singular vectors and values are computed using iterative eigen decomposition of ATA

