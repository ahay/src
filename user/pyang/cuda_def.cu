/* RSF operators using CUDA 
Note: #include "cuda_def.cu" in main function if needed*/
/*
  Copyright (C) 2014 Xi'an Jiaotong University, Pengliang Yang
  
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

#include <cuda_runtime.h>

extern "C" {
#include <rsf.h>
}

#ifndef MAX
#define	MAX(x,y) ((x) > (y) ? (x) : (y))
#endif
#ifndef MIN
#define	MIN(x,y) ((x) < (y) ? (x) : (y))
#endif
#ifndef true
#define true    (1)
#endif
#ifndef false
#define false   (0)
#endif
#ifndef EPS
#define EPS	SF_EPS
#endif
#define Block_Size 512

#define div_up(a, b) (((a) + (b) - 1) / (b))

void sf_check_gpu_error (const char *msg) 
/*< check GPU errors >*/
{
    cudaError_t err = cudaGetLastError ();
    if (cudaSuccess != err) { 
	sf_error ("Cuda error: %s: %s", msg, cudaGetErrorString (err)); 
	exit(0);   
    }
}

/* kernels of basic linear algebric subprograms (BLAS) with CUDA */
__global__ void cuda_scale(float *x, float *y, float alpha, int n)
/*< scale vector: y[]=alpha*x[] >*/
{
	int id=threadIdx.x+blockDim.x*blockIdx.x;
    	if (id<n) y[id]=alpha*x[id];
}

__global__ void cuda_axpy(float *x, float *y, float alpha, int n)
/*< y[]=alpha*x[]+y[] >*/
{
	int id=threadIdx.x+blockDim.x*blockIdx.x;
    	if (id<n) y[id]+=alpha*x[id];
}

__global__  void cuda_dot(float *x, float *y, int n, float *result)
/*< dot product <x,y>:  <<<1, Block_Size>>> >*/
{
  	__shared__ float sdata[Block_Size];
    	int tid=threadIdx.x;
    	sdata[tid]=0.0f;
	for(int s=0; s<(n+Block_Size-1)/Block_Size; s++)
	{
		int id=s*blockDim.x+threadIdx.x;
		sdata[tid] += (id<n)?x[id]*y[id]:0.0;	
	} 
    	__syncthreads();

    	/* do reduction in shared mem */
    	for(int s=blockDim.x/2; s>32; s>>=1) 
    	{
		if (threadIdx.x < s) sdata[tid] += sdata[tid + s]; 
		__syncthreads();
    	}
   	if (tid < 32)
   	{
		if (blockDim.x >=  64) { sdata[tid] += sdata[tid + 32]; }
		if (blockDim.x >=  32) { sdata[tid] += sdata[tid + 16]; }
		if (blockDim.x >=  16) { sdata[tid] += sdata[tid +  8]; }
		if (blockDim.x >=   8) { sdata[tid] += sdata[tid +  4]; }
		if (blockDim.x >=   4) { sdata[tid] += sdata[tid +  2]; }
		if (blockDim.x >=   2) { sdata[tid] += sdata[tid +  1]; }
    	}
     
    	if (tid == 0) { *result=sdata[0]; }
}

__global__ void cuda_mul(bool add, float *z, float *x, float *y, int n)
/*< vector multiplication >*/
{
	int id=threadIdx.x+blockDim.x*blockIdx.x;
    	if (id<n) z[id]=x[id]*y[id]+(add)?z[id]:0.0;
}

