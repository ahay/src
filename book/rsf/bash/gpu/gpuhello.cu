/* GPU example. */
/*
  Copyright (C) 2010 University of Texas at Austin
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include <math.h>
#include <cuda.h>
#include <cuda_runtime_api.h>
extern "C" {
#include <rsf.h>
}
/*#include "kernel.cu"*/

static void sf_check_gpu_error (const char *msg) {
    cudaError_t err = cudaGetLastError ();
    if (cudaSuccess != err) 
        sf_error ("Cuda error: %s: %s", msg,
                   cudaGetErrorString (err));
}

#define BLOCK_SIZE 128
__global__ void gpu_vec_sum (float *a, float *b,
                             float *c) {
    const unsigned int j = blockIdx.x*blockDim.x +
                           threadIdx.x;
    c[j] = a[j] + b[j];
}
int main (int argc, char* argv[]) {
    int n1, n2, esize, i;
    float *a, *b, *c;
    sf_file ain, bin, cout = NULL;
    dim3 dimgrid (1, 1, 1); /* GPU grid */
    dim3 dimblock (BLOCK_SIZE, 1, 1); /* GPU block */
    float *d_a, *d_b, *d_c; /* GPU pointers */

    sf_init (argc, argv);
    cuInit (0); /* Use first GPU device */
    sf_check_gpu_error ("Device initialization");
    cudaSetDevice (0);
    ain = sf_input ("in"); /* Input vector a */
    bin = sf_input ("b"); /* Input vector b */
    if (SF_FLOAT != sf_gettype (ain) ||
        SF_FLOAT != sf_gettype (bin))
        sf_error ("Need float");
    /* Size of an element */
    if (!sf_histint (ain, "esize", &esize))
        esize = sizeof(float);
    /* Vector size */
    if (!sf_histint (ain, "n1", &n1)) sf_error ("No n1=");
    /* Number of vectors */
    n2 = sf_leftsize (ain, 1);
    /* Output vector */
    cout = sf_output ("out");
    /* Vectors in CPU memory */
    a = sf_floatalloc (n1); b = sf_floatalloc (n1);
    c = sf_floatalloc (n1);
    /* Vectors in GPU memory */
    cudaMalloc ((void**)&d_a, n1*esize);
    cudaMalloc ((void**)&d_b, n1*esize);
    cudaMalloc ((void**)&d_c, n1*esize);
    sf_check_gpu_error ("GPU mallocs");
    /* Kernel configuration for this data */
    dimgrid = dim3 (n1/BLOCK_SIZE, 1, 1);
    /* Outer loop over vectors */
    for (i = 0; i < n2; i++) {
        sf_floatread (a, n1, ain); /* Input */
        sf_floatread (b, n1, bin);
        cudaMemcpy (d_a, a, n1*esize, /* a -> GPU */
                    cudaMemcpyHostToDevice);
        cudaMemcpy (d_b, b, n1*esize, /* b -> GPU */
                    cudaMemcpyHostToDevice);
        sf_check_gpu_error ("Copying a&b to GPU");
        /* Parallel summation on GPU */
        gpu_vec_sum<<<dimgrid, dimblock>>>(d_a, d_b, d_c);
        sf_check_gpu_error ("Kernel execution");
        cudaMemcpy (c, d_c, n1*esize, /* GPU -> c */
                    cudaMemcpyDeviceToHost);
        sf_check_gpu_error ("Copying c from GPU");
        sf_floatwrite (c, n1, cout); /* Output */
    }
    sf_fileclose (ain);
    sf_fileclose (bin);
    sf_fileclose (cout);
    return 0;
}

