/* Prestack time migration (2-D/3-D) CUDA kernel. */
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

#ifndef _GPU_KTMIG_KERNEL_H_
#define _GPU_KTMIG_KERNEL_H_

#define BLOCK_SIZE 128
#define MAX_TEXLEN 1 << 27
#define MAX_TEXSIDE 32768
#define MAX_GRID_SIZE 65535

__constant__ int c_obn; /* Number of input traces per each kernel call */
__constant__ int c_nt; /* Number of samples in each input trace */
__constant__ float c_ot; /* Input data t0 */
__constant__ float c_dt; /* Input data sample rate */
__constant__ float c_idt; /* Inverse of input data sample rate */
__constant__ int c_ont; /* Number of samples in each output trace */
__constant__ float c_oot; /* Output data t0 */
__constant__ float c_odt; /* Output data sample rate */

__constant__ float c_oox; /* First x in output */
__constant__ float c_odx; /* Output sample rate in x */
__constant__ int c_onx; /* Number of samples x in output */
__constant__ float c_ooy; /* First y in output */
__constant__ float c_ody; /* Output sample rate in y */
__constant__ int c_ony; /* Number of samples y in output */

__constant__ int c_aox; /* Apperture first indices in x */
__constant__ int c_aoy; /* Apperture first indices in y */
__constant__ int c_anx; /* Apperture number of indices in x */
__constant__ int c_any; /* Apperture number of indices in y */

__constant__ int c_ibn; /* Number of calculated distances per trace */

__constant__ float c_trf; /* Trace factor for antialiasing */
__constant__ float c_trm; /* Maximum half-length of the filter */
__constant__ float c_maxnt; /* Maximum input sample index, usually nt - 1 */

/* Array of shot coordinates for input traces */
texture<float2, 1, cudaReadModeElementType> t_sxy;
/* Array of receiver coordinates for input traces */
texture<float2, 1, cudaReadModeElementType> t_gxy;
/* Array of surface vectors to source and receiver */
texture<float4, 1, cudaReadModeElementType> t_ixy;
/* Input data vector as a texture */
texture<float, 2, cudaReadModeElementType> t_i;
/* Aperture indices */
texture<uint, 1, cudaReadModeElementType> t_ap;

/****************************************************************
 *
 * Differentiation (backward second order) kernel
 *
 ****************************************************************/

/*
  t[gridDim.x*c_nt] - output vector of differentiated traces,
                      also input
*/
__global__ void sf_gpu_ktmig_sbdiff (float *t, const unsigned int n) {
    __shared__ float buffer[BLOCK_SIZE + 2];
    unsigned int tidx = blockIdx.x*c_nt + threadIdx.x;
    float val0 = 0, val1 = 0, val2 = 0;
    int i;

    buffer[threadIdx.x] = 0;
    __syncthreads ();
    for (i = 0; i < n; i++) {
        val0 = t[tidx];
        /* Value at t0 */
        buffer[2 + threadIdx.x] = val0;
        __syncthreads ();
        /* Value at t-1 */
        val1 = buffer[1 + threadIdx.x];
        /* Value at t-2 */
        val2 = buffer[threadIdx.x];
        /* Derivative */
        t[tidx] = 0.5f*c_idt*(3.0f*val0 - 4.0f*val1 + val2);
        __syncthreads ();
        /* Shift everything down */
        buffer[(threadIdx.x + 2)%blockDim.x] = val0;
        __syncthreads ();
        tidx += blockDim.x;
    }
}

/****************************************************************
 *
 * Causal integration kernel
 *
 ****************************************************************/

/*
  t[gridDim.x*c_nt] - output vector of causally integrated traces,
                      also input
*/
__global__ void sf_gpu_ktmig_cint (float *t, const unsigned int n) {
    __shared__ float buffer[BLOCK_SIZE*2];
    unsigned int tidx = blockIdx.x*c_nt + threadIdx.x;
    float val1 = 0, val2 = 0;
    int i, j;

    val1 = t[tidx];
    buffer[threadIdx.x] = val1;
    for (i = 0; i < n; i++) {
        /* Read next portion to the buffer */
        if (i != (n - 1)) {
            val2 = t[tidx + blockDim.x];
            buffer[blockDim.x + threadIdx.x] = val2;
            __syncthreads ();
        }
        /* Integrate causally the current portion */
        for (j = 1; j <= blockDim.x; j++) {
            buffer[threadIdx.x + j] += val1;
            __syncthreads ();
        }
        t[tidx] = buffer[threadIdx.x];
        /* Shift buffer down */
        if (i != (n - 1)) {
            buffer[threadIdx.x] = buffer[blockDim.x + threadIdx.x];
            __syncthreads ();
        }
        tidx += blockDim.x;
        val1 += val2;
    }
}

/****************************************************************
 *
 * Anti-causal integration kernel
 *
 ****************************************************************/

/*
  t[gridDim.x*c_nt] - output vector of anti-causally integrated traces,
                      also input
*/
__global__ void sf_gpu_ktmig_acint (float *t, const unsigned int n) {
    __shared__ float buffer[BLOCK_SIZE*2];
    /* Base index in the output vector */
    unsigned int tidx = (blockIdx.x + 1)*c_nt - blockDim.x + threadIdx.x;
    /* Base index in the local buffer */
    unsigned int ltidx = blockDim.x - threadIdx.x - 1;
    float val1 = 0, val2 = 0;
    int i, j;

    val1 = t[tidx];
    buffer[ltidx] = val1;
    for (i = 0; i < n; i++) {
        /* Read next portion to the buffer */
        if (i != (n - 1)) {
            val2 = t[tidx - blockDim.x];
            buffer[blockDim.x + ltidx] = val2;
            __syncthreads ();
        }
        /* Integrate anti-causally the current portion */
        for (j = 1; j <= blockDim.x; j++) {
            buffer[ltidx + j] += val1;
            __syncthreads ();
        }
        t[tidx] = buffer[ltidx];
        /* Shift buffer down */
        if (i != (n - 1)) {
            buffer[ltidx] = buffer[blockDim.x + ltidx];
            __syncthreads ();
        }
        tidx -= blockDim.x;
        val1 += val2;
    }
}

/****************************************************************
 *
 * Source/receiver-to-image point distance calculation kernel
 *
 ****************************************************************/

/*
  ixy[onx*ony*c_obn] - surface vectors to source and receiver for each image location
*/
__global__ void sf_gpu_ktmig_ixy (float4 *ixy) {
    int i;
    float x, y; /* Image locations */
    float2 xy; /* Source/receiver locations */
    float4 dist;

    /* Image coordinates */
    x = c_oox + (c_aox + (blockIdx.x*blockDim.x + threadIdx.x)%c_anx)*c_odx;
    y = c_ooy + (c_aoy + (blockIdx.x*blockDim.x + threadIdx.x)/c_anx)*c_ody;
    for (i = 0; i < c_obn; i++) {
        /* Source surface vector components */
        xy = tex1Dfetch (t_sxy, i);
        dist.x = x - xy.x;
        dist.y = y - xy.y;
        /* Receiver surface vector components */
        xy = tex1Dfetch (t_gxy, i);
        dist.z = x - xy.x;
        dist.w = y - xy.y;
        ixy[(i*gridDim.x + blockIdx.x)*blockDim.x + threadIdx.x] = dist;
        __syncthreads ();
    }
}

/****************************************************************
 *
 * PSTM kernel with anti-aliasing after Lumley-Claerbout
 *
 ****************************************************************/

/*
  vrms[onx*ony*ont]  - RMS velocity vector,
  image[onx*ony*ont] - output image vector,
  lshift             - shift from the beginning of the array of aperture indices.
*/
__global__ void sf_gpu_ktmig_kernel (float *vrms, float *image,
                                     const unsigned int lshift) {
    /* Index within the output trace */
    const unsigned int tidx = blockDim.x*blockIdx.x + threadIdx.x;
    /* Index within the aperture block */
    const unsigned int bidx = tex1Dfetch (t_ap, blockIdx.y + lshift);
    /* Index within the output vector:
       (Position in y within the grid + position in x within the grid)
       *(number of time samples per output trace)
       + shift in trace */
    const unsigned int oidx = ((c_aoy + bidx/c_anx)*c_onx + c_aox + bidx%c_anx)*c_ont + tidx;
    /* RMS velocity at image location */
    const float v = vrms[oidx];
    /* Slowness at image location */
    const float inv = 1.0f/v;
    const float inv2trf = c_trf*inv*inv;
    /* Pseudodepth^2 at image location */
//  const float depth2 = powf (0.5f*v*(c_oot + tidx*c_odt), 2.0f);
    const float depth2 = 0.25f*v*v*(c_oot + tidx*c_odt)*(c_oot + tidx*c_odt);
    int i;
    float j, k;
    float img = 0.0f, scale, smp;
    float2 vec1;
    float4 vec2;

    /* Loop over input traces */
    for (i = 0; i < c_obn; i++) {
        vec2 = tex1Dfetch (t_ixy, i*c_ibn + bidx);
        /* vec1.x - squared distance to source from the image point on the surface,
           vec1.y - squared distance to receiver from the image point on the surface */
        vec1.x = vec2.x*vec2.x + vec2.y*vec2.y;
        vec1.y = vec2.z*vec2.z + vec2.w*vec2.w;
        /* Time from source to image point in pseudodepth */
        vec1.x = sqrtf (vec1.x + depth2)*inv;
        /* Time from receiver to image point in pseudodepth */
        vec1.y = sqrtf (vec1.y + depth2)*inv;
        /* double root square time = time to source + time to receiver */
        j = (vec1.x + vec1.y - c_ot)*c_idt; /* Input sample index */
        /* (distance to source.x)/(time to source) + (distance to receiver.x)/(time to receiver) */
        vec2.x = vec2.x/vec1.x + vec2.z/vec1.y;
        /* (distance to source.y)/(time to source) + (distance to receiver.y)/(time to receiver) */
        vec2.y = vec2.y/vec1.x + vec2.w/vec1.y;
        /* Filter length */
        k = inv2trf*sqrtf (vec2.x*vec2.x + vec2.y*vec2.y);
        /* Truncate filter */
        k = fminf (k, c_trm);
        /* If any of the three points is out of range - zero everything out */
        if ((j - k - 1.0f) >= 0.0f && (j + k + 1.0f) <= c_maxnt) {
            /* Scaling factor */
            scale = 1.0f/(1.0f + k);
            scale *= scale;
            /* Collect samples */
            smp = 2.0f*tex2D (t_i, j + 0.5f, i + 0.5f)
                  -tex2D (t_i, j - k - 0.5f, i + 0.5f)
                  -tex2D (t_i, j + k + 1.5f, i + 0.5f);
            /* Contribute to the image point */
            img += scale*smp;
        }
    }

    image[oidx] += img;
}

/****************************************************************
 *
 * PSTM kernel without anti-aliasing (simplified version from above)
 *
 ****************************************************************/

/*
  vrms[onx*ony*ont]  - RMS velocity vector,
  image[onx*ony*ont] - output image vector,
  lshift             - shift from the beginning of the array of aperture indices
*/
__global__ void sf_gpu_ktmig_noaa_kernel (float *vrms, float *image,
                                          const unsigned int lshift) {
    /* Index within the output trace */
    const unsigned int tidx = blockDim.x*blockIdx.x + threadIdx.x;
    /* Index within the aperture block */
    const unsigned int bidx = tex1Dfetch (t_ap, blockIdx.y + lshift);
    /* Index within the output vector:
       (Position in y within the grid + position in x within the grid)
       *(number of time samples per output trace)
       + shift in trace */
    const unsigned int oidx = ((c_aoy + bidx/c_anx)*c_onx + c_aox + bidx%c_anx)*c_ont + tidx;
    /* RMS velocity at image location */
    const float v = vrms[oidx];
    /* Slowness at image location */
    const float inv = 1.0f/v;
    /* Pseudodepth^2 at image location */
//  const float depth2 = powf (0.5f*v*(c_oot + tidx*c_odt), 2.0f);
    const float depth2 = 0.25f*v*v*(c_oot + tidx*c_odt)*(c_oot + tidx*c_odt);
    int i;
    float j, k;
    float img = 0.0f;
    float2 vec1;
    float4 vec2;

    /* Loop over input traces */
    for (i = 0; i < c_obn; i++) {
        vec2 = tex1Dfetch (t_ixy, i*c_ibn + bidx);
        /* vec1.x - squared distance to source from the image point on the surface,
           vec1.y - squared distance to receiver from the image point on the surface */
        vec1.x = vec2.x*vec2.x + vec2.y*vec2.y;
        vec1.y = vec2.z*vec2.z + vec2.w*vec2.w;
        /* Time from source to image point in pseudodepth */
        vec1.x = sqrtf (vec1.x + depth2)*inv;
        /* Time from receiver to image point in pseudodepth */
        vec1.y = sqrtf (vec1.y + depth2)*inv;
        /* double root square time = time to source + time to receiver */
        j = (vec1.x + vec1.y - c_ot)*c_idt; /* Input sample index */
        if (j >= 0.0f && j <= c_maxnt) {
            k = (float)i + 0.5f;
            /* Contribute to the image point */
            img += tex2D (t_i, j + 0.5f, k);
        }
    }

    image[oidx] += img;
}

#endif /* _GPU_KTMIG_KERNEL_H_ */

