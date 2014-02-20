/* 2D prestack GPU-based RTM using effective boundary saving
Some basic descriptions of this code are in order.
1) Coordinate configuration of seismic data:

		o--------------> x (2nd dim: *.y)
		|
		|
		|
		|
		|
		z (1st dim: *.x)
	1st dim: i1=threadIdx.x+blockDim.x*blockIdx.x;
	2nd dim: i2=threadIdx.y+blockDim.y*blockIdx.y;
	(i1, i2)=i1+i2*nnz;

2) stability condition:	
	min(dx, dz)>sqrt(2)*dt*max(v) (NJ=2)
   numerical dispersion condition:	
	max(dx, dz)<min(v)/(10*fmax)  (NJ=2)
	max(dx, dz)<min(v)/(5*fmax)   (NJ=4)

3) This code doesn't save the history of forward time steps. We 
   just save the least boundaries (referred to as effective boundary 
   in our work) of every time step and the two final steps of the 
   wavefield. Using this information, we can easily reconstruct 
   the exact wavefield in the reverse time steps. It is noteworthy
   that to implement large scale seismic imaging, pinned memory is 
   employed to save the boundaries of each step so that all the saved
   data can be computed on the device directly.

4) In our implementation, we employ staggered grid based 
   convolutional PML (CPML) boundary condition. Using 20 points for 
   CPML is enough to obtain perfect absorbing effect (while commonly 
   used sponge ABC may need 30 or more). However, we use 32 points on
   each side due to the grid alignment reasons. (To make your code 
   fast, you should consider that the GPU codes implementation unit 
   is half-warp (16 threads). The thickness of the boundary should be 
   times of 16. 

5) The final images can be two kinds: result of correlation imaging 
   condition and the normalized one. The normalized correlation imaging
   result is preferred due to compensated illumination. This code does
   not perform any kind of filtering, which is recommended if you 
   obtained the CUDA RTM result. Some of the filters are popular and 
   effective to remove the low frequency artifacts of the imaging: the
   Laplacian filtering, derivative filtering and the bandpass filtering. 
   In this code, we use laplacian filtering.
*/

/*
  Copyright (C) 2013  Xi'an Jiaotong University (Pengliang Yang)
    Email: ypl.2100@gmail.com	    
  Acknowledgement: This code is written with the help of Baoli Wang.

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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
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
#define EPS	1.0e-15f
#endif
#ifndef PI
#define PI 	3.141592653589793f
#endif
#define Block_Size1 16		// 1st dim block size
#define Block_Size2 16		// 2nd dim block size
const int npml=32;		// thickness of PML boundary
const int nbell=3;		// radius of Gaussian bell 
//const bool frsf=fasle;	// free surface on top boundary or not

#include "cuda_kernels.cu"

static bool 	csdgather; 	// common shot gather (CSD) or not 
static int 	nz1,nx1, nz, nx, nnz, nnx, N, NJ, ns, ng, nt, nt_h;
static int 	jsx,jsz,jgx,jgz,sxbeg,szbeg,gxbeg,gzbeg;
static float 	fm, dt, dz, dx, _dz, _dx, vmute;
static dim3 	dimbbell, dimg0, dimb0;
static dim3 	dimglr1, dimblr1, dimglr2, dimblr2;//lr=left and right
static dim3 	dimgtb1, dimbtb1, dimgtb2, dimbtb2;//tb=top and bottom

// variables on host
float 	*seis, *v0, *vel, *p;
// variables on device
int 	*d_Sxz, *d_Gxz;				// set source and geophone position
float 	*d_bell,*d_wlt, *d_dobs,  *d_vel;	// bell, wavelet, seismograms, velocity (vel)
float 	*d_sp0, *d_sp1, *d_sp2, *d_svx, *d_svz;	// p, vx, vz for sources
float 	*d_gp0, *d_gp1, *d_gp2, *d_gvx, *d_gvz;	// p, vx, vz for geophones
float 	*d_bx1, *d_bx2, *d_bz1, *d_bz2;		// PML ABC coefficients for p and v (vx, vz)
float 	*d_convpx, *d_convpz, *d_convvx, *d_convvz;// auxiliary variables to decay p and v in PML zone
float 	*d_Iss, *d_Isg, *d_I1,*d_I2;		// I1: image without normalization; I2: normalized image; 
float 	*h_boundary, *d_boundary;		// boundary on host and device
float	*ptr=NULL;

void matrix_transpose(float *matrix, int nx, int nz)
{
	float *tmp=(float*)malloc(nx*nz*sizeof(float));
	if (tmp==NULL) {sf_warning("out of memory!"); exit(1);}
	for(int iz=0; iz<nz; iz++){
		for(int ix=0; ix<nx; ix++){
			tmp[iz+nz*ix]=matrix[ix+nx*iz];
		}
	}
	memcpy(matrix, tmp, nx*nz*sizeof(float));
	free(tmp);
}



void expand(float *b, float *a, int npml, int nnz, int nnx, int nz1, int nx1)
/*< expand domain of 'a' to 'b':  a, size=nz1*nx1; b, size=nnz*nnx;  >*/
{
    int iz,ix;
    for     (ix=0;ix<nx1;ix++) {
	for (iz=0;iz<nz1;iz++) {
	    b[(npml+ix)*nnz+(npml+iz)] = a[ix*nz1+iz];
	}
    }
    for     (ix=0; ix<nnx; ix++) {
	for (iz=0; iz<npml; iz++)   	b[ix*nnz+iz] = b[ix*nnz+npml];//top
	for (iz=nz1+npml; iz<nnz; iz++) b[ix*nnz+iz] = b[ix*nnz+npml+nz1-1];//bottom
    }

    for (iz=0; iz<nnz; iz++){
	for(ix=0; ix<npml; ix++) 	b[ix*nnz+iz] = b[npml*nnz+iz];//left
	for(ix=npml+nx1; ix<nnx; ix++)	b[ix*nnz+iz] = b[(npml+nx1-1)*nnz+iz];//right
    }
}

void window(float *to, float *from, int npml, int nnz, int nnx, int nz1, int nx1)
/*< window from: size=nnz*nnx; to: size=nz1*nx1; >*/
{
	int ix, iz;
	for(ix=0;ix<nx1;ix++){
		for(iz=0;iz<nz1;iz++){
			to[iz+ix*nz1]=from[(iz+npml)+(ix+npml)*nnz];
		}
	}
}


void check_grid_sanity(int NJ, float *vel, float fm, float dz, float dx, float dt, int N)
/*< sanity check about stability condition and non-dispersion condition >*/
{
	float C;
	if(NJ==2) C=1;
	else if (NJ==4)	 	C=0.857;
	else if (NJ==6)		C=0.8;
	else if (NJ==8) 	C=0.777;
	else if (NJ==10)	C=0.759;

	float maxvel=vel[0], minvel=vel[0];
	for(int i=0; i<N; i++)	{
		if(vel[i]>maxvel) maxvel=vel[i];
		if(vel[i]<minvel) minvel=vel[i];
	}
	float tmp=dt*maxvel*sqrtf(1.0/(dx*dx)+1.0/(dz*dz));

	if (tmp>=C) sf_warning("Stability condition not satisfied!");
	if (fm>=minvel/(5*MAX(dx,dz))) sf_warning("Non-dispersion relation not satisfied!");
}


void device_alloc()
/*< allocate variable space on device >*/
{
	cudaMalloc(&d_bell,	(2*nbell+1)*(2*nbell+1)*sizeof(float));
    	cudaMalloc(&d_Sxz,	ns*sizeof(int));
    	cudaMalloc(&d_Gxz, 	ng*sizeof(int));
	cudaMalloc(&d_wlt,	nt*sizeof(float));
    	cudaMalloc(&d_dobs, 	ng*nt*sizeof(float));
    	cudaMalloc(&d_vel, 	N*sizeof(float));
    	cudaMalloc(&d_sp0, 	N*sizeof(float));
    	cudaMalloc(&d_sp1, 	N*sizeof(float));
    	cudaMalloc(&d_sp2, 	N*sizeof(float));
    	cudaMalloc(&d_svx, 	N*sizeof(float));
    	cudaMalloc(&d_svz, 	N*sizeof(float));
    	cudaMalloc(&d_gp0, 	N*sizeof(float));
    	cudaMalloc(&d_gp1, 	N*sizeof(float));
    	cudaMalloc(&d_gp2, 	N*sizeof(float));
    	cudaMalloc(&d_gvx, 	N*sizeof(float));
    	cudaMalloc(&d_gvz, 	N*sizeof(float));
    	cudaMalloc(&d_bx1, 	2*npml*nnz*sizeof(float));// left and right ABC coefficients for p
    	cudaMalloc(&d_bz1, 	2*npml*nnx*sizeof(float));// top and bottom ABC coefficients for p
    	cudaMalloc(&d_bx2, 	2*npml*nnz*sizeof(float));// left and right ABC coefficients for v (vx, vz)
    	cudaMalloc(&d_bz2, 	2*npml*nnx*sizeof(float));// top and bottom ABC coefficients for v (vx, vz)
    	cudaMalloc(&d_convpx, 	2*npml*nnz*sizeof(float));//  (left and right)
    	cudaMalloc(&d_convpz, 	2*npml*nnx*sizeof(float));//  (top and bottom)
    	cudaMalloc(&d_convvx, 	2*npml*nnz*sizeof(float));//  (left and right)
    	cudaMalloc(&d_convvz, 	2*npml*nnx*sizeof(float));//  (top and bottom)
    	cudaMalloc(&d_Iss, 	N*sizeof(float));
    	cudaMalloc(&d_Isg, 	N*sizeof(float));
    	cudaMalloc(&d_I1, 	N*sizeof(float));
    	cudaMalloc(&d_I2, 	N*sizeof(float));
	cudaHostAlloc(&h_boundary, nt_h*2*(NJ-1)*(nx+nz)*sizeof(float), cudaHostAllocMapped);	
	cudaMalloc(&d_boundary, (nt-nt_h)*2*(NJ-1)*(nx+nz)*sizeof(float));


    	cudaError_t err = cudaGetLastError ();
    	if (cudaSuccess != err) 
	sf_warning("Cuda error: Failed to allocate required memory!: %s", cudaGetErrorString(err));
}


void device_free()
/*< free the variables oon device >*/
{
	cudaFree(d_bell);
    	cudaFree(d_Sxz);
    	cudaFree(d_Gxz);
	cudaFree(d_wlt);
    	cudaFree(d_dobs);
    	cudaFree(d_vel);
    	cudaFree(d_sp0);
    	cudaFree(d_sp1);
    	cudaFree(d_sp2);
    	cudaFree(d_svx);
    	cudaFree(d_svz);
    	cudaFree(d_gp0);
    	cudaFree(d_gp1);
    	cudaFree(d_gp2);
    	cudaFree(d_gvx);
    	cudaFree(d_gvz);
	cudaFree(d_bx1);
	cudaFree(d_bx2);
	cudaFree(d_bz1);
	cudaFree(d_bz2);
    	cudaFree(d_convpx);
    	cudaFree(d_convpz);
    	cudaFree(d_convvx);
    	cudaFree(d_convvz);
	cudaFree(d_Iss);
	cudaFree(d_Isg);
	cudaFree(d_I1);
	cudaFree(d_I2);
	cudaFreeHost(h_boundary);
    	cudaFree(d_boundary);

    	cudaError_t err = cudaGetLastError ();
    	if (cudaSuccess != err)
	sf_warning("Cuda error: Failed to free the allocated memory!: %s", cudaGetErrorString(err));
}

void wavefield_init(float *d_p0, float *d_p1, float *d_p2, float *d_vx, float *d_vz, float *d_convpx, float *d_convpz, float *d_convvx, float *d_convvz)
/*< initialize wavefield variables >*/
{
	cudaMemset(d_p0, 	0,	N*sizeof(float));
	cudaMemset(d_p1, 	0,	N*sizeof(float));
	cudaMemset(d_p2, 	0,  	N*sizeof(float));
	cudaMemset(d_vx, 	0,	N*sizeof(float));
	cudaMemset(d_vz, 	0,	N*sizeof(float));
	cudaMemset(d_convpx, 	0,	2*npml*nnz*sizeof(float));
	cudaMemset(d_convpz, 	0,	2*npml*nnx*sizeof(float));
	cudaMemset(d_convvx, 	0,	2*npml*nnz*sizeof(float));
	cudaMemset(d_convvz, 	0,	2*npml*nnx*sizeof(float));

    	cudaError_t err = cudaGetLastError ();
    	if (cudaSuccess != err) 
	printf("Cuda error: Failed to initialize the wavefield variables!: %s", cudaGetErrorString(err));
}

void step_forward(float *vel, float *d_p0, float *d_p1, float *d_p2, float *d_vx, float *d_vz, float *d_convvx, float *d_convvz, float *d_convpx, float *d_convpz, float *d_bx1, float *d_bz1, float *d_bx2, float *d_bz2)
/*< step of forward propagation >*/
{
	// p0: p{it-1}; 
	// p1: p{it};
	// p2: p{it+1};
	if (NJ==2)	{
		cuda_forward_v_2<<<dimg0, dimb0>>>(d_p1, d_vx, d_vz, _dx, _dz, npml, nnz, nnx);
		cuda_PML_vz_2<<<dimgtb1, dimbtb1>>>(d_p1, d_convpz, d_bz2, d_vz, _dz, npml, nnz, nnx);
		cuda_PML_vx_2<<<dimglr1, dimblr1>>>(d_p1, d_convpx, d_bx2, d_vx, _dx, npml, nnz, nnx);
		cuda_forward_p_2<<<dimg0, dimb0>>>(d_vel, d_p0, d_p1, d_p2, d_vx, d_vz, dt, _dx, _dz, npml, nnz, nnx);
		cuda_PML_pz_2<<<dimgtb1, dimbtb1>>>(d_vel, d_p2, d_convvz, d_bz1, d_vz, dt, _dz, npml, nnz, nnx);
		cuda_PML_px_2<<<dimglr1, dimblr1>>>(d_vel, d_p2, d_convvx, d_bx1, d_vx, dt, _dx, npml, nnz, nnx);
	}else if (NJ==4){
		cuda_forward_v_4<<<dimg0, dimb0>>>(d_p1, d_vx, d_vz, _dx, _dz, npml, nnz, nnx);
		cuda_PML_vz_4<<<dimgtb1, dimbtb1>>>(d_p1, d_convpz, d_bz2, d_vz, _dz, npml, nnz, nnx);
		cuda_PML_vx_4<<<dimglr1, dimblr1>>>(d_p1, d_convpx, d_bx2, d_vx, _dx, npml, nnz, nnx);
		cuda_forward_p_4<<<dimg0, dimb0>>>(d_vel, d_p0, d_p1, d_p2, d_vx, d_vz, dt, _dx, _dz, npml, nnz, nnx);
		cuda_PML_pz_4<<<dimgtb1, dimbtb1>>>(d_vel, d_p2, d_convvz, d_bz1, d_vz, dt, _dz, npml, nnz, nnx);
		cuda_PML_px_4<<<dimglr1, dimblr1>>>(d_vel, d_p2, d_convvx, d_bx1, d_vx, dt, _dx, npml, nnz, nnx);
	}else if (NJ==6){
		cuda_forward_v_6<<<dimg0, dimb0>>>(d_p1, d_vx, d_vz, _dx, _dz, npml, nnz, nnx);
		cuda_PML_vz_6<<<dimgtb1, dimbtb1>>>(d_p1, d_convpz, d_bz2, d_vz, _dz, npml, nnz, nnx);
		cuda_PML_vx_6<<<dimglr1, dimblr1>>>(d_p1, d_convpx, d_bx2, d_vx, _dx, npml, nnz, nnx);
		cuda_forward_p_6<<<dimg0, dimb0>>>(d_vel, d_p0, d_p1, d_p2, d_vx, d_vz, dt, _dx, _dz, npml, nnz, nnx);
		cuda_PML_pz_6<<<dimgtb1, dimbtb1>>>(d_vel, d_p2, d_convvz, d_bz1, d_vz, dt, _dz, npml, nnz, nnx);
		cuda_PML_px_6<<<dimglr1, dimblr1>>>(d_vel, d_p2, d_convvx, d_bx1, d_vx, dt, _dx, npml, nnz, nnx);
	}else if (NJ==8){
		cuda_forward_v_8<<<dimg0, dimb0>>>(d_p1, d_vx, d_vz, _dx, _dz, npml, nnz, nnx);
		cuda_PML_vz_8<<<dimgtb1, dimbtb1>>>(d_p1, d_convpz, d_bz2, d_vz, _dz, npml, nnz, nnx);
		cuda_PML_vx_8<<<dimglr1, dimblr1>>>(d_p1, d_convpx, d_bx2, d_vx, _dx, npml, nnz, nnx);
		cuda_forward_p_8<<<dimg0, dimb0>>>(d_vel, d_p0, d_p1, d_p2, d_vx, d_vz, dt, _dx, _dz, npml, nnz, nnx);
		cuda_PML_pz_8<<<dimgtb1, dimbtb1>>>(d_vel, d_p2, d_convvz, d_bz1, d_vz, dt, _dz, npml, nnz, nnx);
		cuda_PML_px_8<<<dimglr1, dimblr1>>>(d_vel, d_p2, d_convvx, d_bx1, d_vx, dt, _dx, npml, nnz, nnx);
	}else if (NJ==10){
		cuda_forward_v_10<<<dimg0, dimb0>>>(d_p1, d_vx, d_vz, _dx, _dz, npml, nnz, nnx);
		cuda_PML_vz_10<<<dimgtb1, dimbtb1>>>(d_p1, d_convpz, d_bz2, d_vz, _dz, npml, nnz, nnx);
		cuda_PML_vx_10<<<dimglr1, dimblr1>>>(d_p1, d_convpx, d_bx2, d_vx, _dx, npml, nnz, nnx);
		cuda_forward_p_10<<<dimg0, dimb0>>>(d_vel, d_p0, d_p1, d_p2, d_vx, d_vz, dt, _dx, _dz, npml, nnz, nnx);
		cuda_PML_pz_10<<<dimgtb1, dimbtb1>>>(d_vel, d_p2, d_convvz, d_bz1, d_vz, dt, _dz, npml, nnz, nnx);
		cuda_PML_px_10<<<dimglr1, dimblr1>>>(d_vel, d_p2, d_convvx, d_bx1, d_vx, dt, _dx, npml, nnz, nnx);
	}
}

void step_backward(float *d_vel, float *d_p0, float *d_p1, float *d_p2, float *d_vx, float *d_vz)
/*< step of backward propagation >*/
{
	// p0: p{it-1}; 
	// p1: p{it};
	// p2: p{it+1};
	if (NJ==2){
		cuda_forward_v_2<<<dimg0, dimb0>>>(d_p1, d_vx, d_vz, _dx, _dz, npml, nnz, nnx);
		cuda_forward_p_2<<<dimg0, dimb0>>>(d_vel, d_p0, d_p1, d_p2, d_vx, d_vz, dt, _dx, _dz, npml, nnz, nnx);
	}
	else if (NJ==4)	{
		cuda_forward_v_4<<<dimg0, dimb0>>>(d_p1, d_vx, d_vz, _dx, _dz, npml, nnz, nnx);
		cuda_forward_p_4<<<dimg0, dimb0>>>(d_vel, d_p0, d_p1, d_p2, d_vx, d_vz, dt, _dx, _dz, npml, nnz, nnx);
	}
	else if (NJ==6)	{
		cuda_forward_v_6<<<dimg0, dimb0>>>(d_p1, d_vx, d_vz, _dx, _dz, npml, nnz, nnx);
		cuda_forward_p_6<<<dimg0, dimb0>>>(d_vel, d_p0, d_p1, d_p2, d_vx, d_vz, dt, _dx, _dz, npml, nnz, nnx);
	}
	else if (NJ==8)	{
		cuda_forward_v_8<<<dimg0, dimb0>>>(d_p1, d_vx, d_vz, _dx, _dz, npml, nnz, nnx);
		cuda_forward_p_8<<<dimg0, dimb0>>>(d_vel, d_p0, d_p1, d_p2, d_vx, d_vz, dt, _dx, _dz, npml, nnz, nnx);
	}
	else if (NJ==10){
		cuda_forward_v_8<<<dimg0, dimb0>>>(d_p1, d_vx, d_vz, _dx, _dz, npml, nnz, nnx);
		cuda_forward_p_8<<<dimg0, dimb0>>>(d_vel, d_p0, d_p1, d_p2, d_vx, d_vz, dt, _dx, _dz, npml, nnz, nnx);
	}
}


int main(int argc, char* argv[])
{
	float phost;
	int tdmute;
    	sf_file vmodl, imag1, imag2; /* I/O files */

    	/* initialize Madagascar */
    	sf_init(argc,argv);

    	/*< set up I/O files >*/
    	vmodl = sf_input ("in");   /* velocity model, unit=m/s */
    	imag1 = sf_output("out");  /* output image with correlation imaging condition */ 
    	imag2 = sf_output("imag2");  /* output image with normalized correlation imaging condition */ 

    	/* get parameters for RTM */
    	if (!sf_histint(vmodl,"n1",&nz1)) sf_error("no n1");
    	if (!sf_histint(vmodl,"n2",&nx1)) sf_error("no n2");
    	if (!sf_histfloat(vmodl,"d1",&dz)) sf_error("no d1");
   	if (!sf_histfloat(vmodl,"d2",&dx)) sf_error("no d2");

    	if (!sf_getfloat("fm",&fm)) sf_error("no fm");	/* dominant freq of ricker */
    	if (!sf_getfloat("dt",&dt)) sf_error("no dt");	/* time interval */

    	if (!sf_getint("nt",&nt))   sf_error("no nt");	/* total modeling time steps */
    	if (!sf_getint("ns",&ns))   sf_error("no ns");	/* total shots */
    	if (!sf_getint("ng",&ng))   sf_error("no ng");	/* total receivers in each shot */
	
    	if (!sf_getint("jsx",&jsx))   sf_error("no jsx");/* source x-axis  jump interval  */
    	if (!sf_getint("jsz",&jsz))   jsz=0;/* source z-axis jump interval  */
    	if (!sf_getint("jgx",&jgx))   jgx=1;/* receiver x-axis jump interval */
    	if (!sf_getint("jgz",&jgz))   jgz=0;/* receiver z-axis jump interval */
    	if (!sf_getint("sxbeg",&sxbeg))   sf_error("no sxbeg");/* x-begining index of sources, starting from 0 */
    	if (!sf_getint("szbeg",&szbeg))   sf_error("no szbeg");/* z-begining index of sources, starting from 0 */
    	if (!sf_getint("gxbeg",&gxbeg))   sf_error("no gxbeg");/* x-begining index of receivers, starting from 0 */
    	if (!sf_getint("gzbeg",&gzbeg))   sf_error("no gzbeg");/* z-begining index of receivers, starting from 0 */

    	if (!sf_getint("order",&NJ))   NJ=6;/* order of finite difference, order=2,4,6,8,10 */
    	if (!sf_getfloat("phost",&phost)) phost=0;/* phost% points on host with zero-copy pinned memory, the rest on device */
	if (!sf_getbool("csdgather",&csdgather)) csdgather=true;/* default, common shot-gather; if n, record at every point*/
	if (!sf_getfloat("vmute",&vmute))   vmute=1500;/* muting velocity to remove the low-freq artifacts, unit=m/s*/
	if (!sf_getint("tdmute",&tdmute))   tdmute=200;/* number of deleyed time samples to mute */

    	sf_putint(imag1,"n1",nz1);
    	sf_putint(imag1,"n2",nx1);
    	sf_putfloat(imag1,"d1",dz);
    	sf_putfloat(imag1,"d2",dx);
    	sf_putint(imag2,"n1",nz1);
    	sf_putint(imag2,"n2",nx1);
    	sf_putfloat(imag2,"d1",dz);
    	sf_putfloat(imag2,"d2",dx);
	
	_dx=1.0/dx;
	_dz=1.0/dz;
	nt_h=0.01*phost*nt; 

	nx=(int)((nx1+Block_Size1-1)/Block_Size1)*Block_Size1;
	nz=(int)((nz1+Block_Size2-1)/Block_Size2)*Block_Size2;
    	nnz = 2*npml+nz;
    	nnx = 2*npml+nx;
	N=nnz*nnx;
	dimbbell=dim3(2*nbell+1,2*nbell+1);
    	dimb0=dim3(Block_Size1, Block_Size2);  	dimg0=dim3(nnz/Block_Size1, nnx/Block_Size2);
	dimblr1=dim3(Block_Size1, 32);		dimglr1=dim3(nnz/Block_Size1, 2);
	dimbtb1=dim3(32, Block_Size2); 		dimgtb1=dim3(2, nnx/Block_Size2);
	dimblr2=dim3(nnz/Block_Size1,(NJ+15)/16); dimglr2=dim3(Block_Size1, 16);
	dimbtb2=dim3(16, Block_Size2);		dimgtb2=dim3((NJ+15)/16, nnx/Block_Size2);
	
    	v0=(float*)malloc(nx1*nz1*sizeof(float));
	if (v0==NULL) 	{ sf_warning("out of memory!"); exit(1);}
    	seis=(float*)malloc(ng*nt*sizeof(float));
	if (seis==NULL) { sf_warning("out of memory!"); exit(1);}
    	vel=(float*)malloc(N*sizeof(float));
	if (vel==NULL) 	{ sf_warning("out of memory!"); exit(1);}
    	p=(float*)malloc(N*sizeof(float));
	if (p==NULL) 	{ sf_warning("out of memory!"); exit(1);}
    	memset(v0, 0, nz1*nx1*sizeof(float));
    	memset(seis, 0, ng*nt*sizeof(float));
    	memset(vel, 0, N*sizeof(float));
    	memset(p, 0, N*sizeof(float));

   	sf_floatread(v0,nz1*nx1,vmodl);
	expand(vel, v0, npml, nnz, nnx, nz1, nx1);
	check_grid_sanity(NJ, vel, fm, dz, dx, dt, N);

    	cudaSetDevice(0);
    	cudaError_t err = cudaGetLastError ();
    	if (cudaSuccess != err) 
	sf_warning("Cuda error: Failed to initialize device: %s", cudaGetErrorString(err));
	device_alloc(); 

	float mstimer = 0;// timer unit: millionseconds
	cudaEvent_t start, stop;
  	cudaEventCreate(&start);	
	cudaEventCreate(&stop);

	cuda_init_bell<<<dim3(1,1),dim3(2*nbell+1,2*nbell+1)>>>(d_bell);
	cuda_ricker_wavelet<<<(nt+511)/512,512>>>(d_wlt, fm, dt, nt);
	if (!(sxbeg>=0 && szbeg>=0 && sxbeg+(ns-1)*jsx<nx && szbeg+(ns-1)*jsz<nz))	
	{ sf_warning("sources exceeds the computing zone!"); exit(1);}
	cuda_set_sg<<<(ns+255)/256, 256>>>(d_Sxz, sxbeg, szbeg, jsx, jsz, ns, npml, nnz);

	int distx=sxbeg-gxbeg;
	int distz=szbeg-gzbeg;
	if (csdgather)	{
		if (!(gxbeg>=0 && gzbeg>=0 && gxbeg+(ng-1)*jgx<nx && gzbeg+(ng-1)*jgz<nz &&
		(sxbeg+(ns-1)*jsx)+(ng-1)*jgx-distx <nx  && (szbeg+(ns-1)*jsz)+(ng-1)*jgz-distz <nz))	
		{ sf_warning("geophones exceeds the computing zone!"); exit(1);}
	}
	else{
		if (!(gxbeg>=0 && gzbeg>=0 && gxbeg+(ng-1)*jgx<nx && gzbeg+(ng-1)*jgz<nz))	
		{ sf_warning("geophones exceeds the computing zone!"); exit(1);}
	}
	cuda_set_sg<<<(ng+255)/256, 256>>>(d_Gxz, gxbeg, gzbeg, jgx, jgz, ng, npml, nnz);

    	cudaMemcpy(d_vel, vel, N*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemset(d_Iss, 0,  	N*sizeof(float));
	cudaMemset(d_Isg, 0,  	N*sizeof(float));
	cudaMemset(d_I1, 0,  	N*sizeof(float));
	cudaMemset(d_I2, 0,  	N*sizeof(float));    	
	cuda_init_abcz<<<dimgtb1, dimbtb1>>>(d_vel, d_bz1, d_bz2, dx, dz, dt, npml, nnz, nnx);
	cuda_init_abcx<<<dimglr1, dimblr1>>>(d_vel, d_bx1, d_bx2, dx, dz, dt, npml, nnz, nnx);

    	for(int is=0; is<ns; is++)
	{
		cudaEventRecord(start);

	    	cudaMemset(d_Isg, 	0,	N*sizeof(float));
		cudaMemset(d_Iss, 	0,	N*sizeof(float));
		cudaMemset(d_dobs, 	0,	nt*ng*sizeof(float));
		cudaMemset(h_boundary, 	0,	nt_h*2*(NJ-1)*(nx+nz)*sizeof(float));
	    	cudaMemset(d_boundary, 	0,	(nt-nt_h)*2*(NJ-1)*(nx+nz)*sizeof(float));
		wavefield_init(d_sp0, d_sp1, d_sp2, d_svx, d_svz, d_convpx, d_convpz, d_convvx, d_convvz);
		if (csdgather)	{
			gxbeg=sxbeg+is*jsx-distx;
			cuda_set_sg<<<(ng+255)/256, 256>>>(d_Gxz, gxbeg, gzbeg, jgx, jgz, ng, npml, nnz);
		}
		for(int kt=0; kt<nt; kt++)
		{
			cuda_add_bellwlt<<<dim3(1,1), dimbbell>>>(d_sp1, d_bell, &d_wlt[kt], &d_Sxz[is], 1, npml, nnz, nnx, true);
			//cuda_add_source<<<1,1>>>(d_sp1, &d_wlt[kt], &d_Sxz[is], 1, true);
			step_forward(d_vel, d_sp0, d_sp1, d_sp2, d_svx, d_svz, d_convvx, d_convvz, d_convpx, d_convpz, d_bx1, d_bz1, d_bx2, d_bz2);
			ptr=d_sp0; d_sp0=d_sp1; d_sp1=d_sp2; d_sp2=ptr;

			cuda_record<<<(ng+255)/256, 256>>>(d_sp0, &d_dobs[kt*ng], d_Gxz, ng);
			cuda_mute<<<(ng+511)/512, 512>>>(&d_dobs[kt*ng], gzbeg, szbeg, gxbeg, sxbeg+is*jsx, jgx, kt, tdmute, vmute, dt, dz, dx, ng);

			if(kt<nt_h) cudaHostGetDevicePointer(&ptr, &h_boundary[kt*2*(NJ-1)*(nx+nz)], 0);
			else  ptr=&d_boundary[(kt-nt_h)*2*(NJ-1)*(nx+nz)];
			cuda_rw_innertb<<<dimgtb2, dimbtb2>>>(ptr, 	d_sp0, npml, nnz, nnx, NJ, false);
			cuda_rw_innerlr<<<dimglr2, dimblr2>>>(&ptr[2*(NJ-1)*nx], d_sp0, npml, nnz, nnx, NJ, false);
		}

		ptr=d_sp0; d_sp0=d_sp1; d_sp1=ptr;
		wavefield_init(d_gp0, d_gp1, d_gp2, d_gvx, d_gvz, d_convpx, d_convpz, d_convvx, d_convvz);
		for(int kt=nt-1; kt>-1; kt--)
		{
			// read saved boundary
			if(kt<nt_h) cudaHostGetDevicePointer(&ptr, &h_boundary[kt*2*(NJ-1)*(nx+nz)], 0);
			else  ptr=&d_boundary[(kt-nt_h)*2*(NJ-1)*(nx+nz)];
			cuda_rw_innertb<<<dimgtb2, dimbtb2>>>(ptr, 		d_sp1, npml, nnz, nnx, NJ, true);
			cuda_rw_innerlr<<<dimglr2, dimblr2>>>(&ptr[2*(NJ-1)*nx], d_sp1, npml, nnz, nnx, NJ, true);

			// backward time step source wavefield
			step_backward(d_vel, d_sp0, d_sp1, d_sp2, d_svx, d_svz);
			// subtract the wavelet
			cuda_add_bellwlt<<<dim3(1,1), dimbbell>>>(d_sp1, d_bell, &d_wlt[kt], &d_Sxz[is], 1, npml, nnz, nnx, false);
			//cuda_add_source<<<1,1>>>(d_sp1, &d_wlt[kt], &d_Sxz[is], 1, false);
			ptr=d_sp0; d_sp0=d_sp1; d_sp1=d_sp2; d_sp2=ptr;

			// add receiver term
			cuda_add_source<<<(ng+255)/256,256>>>(d_gp0, &d_dobs[kt*ng], d_Gxz, ng, true);
			// backward time step receiver wavefield
			step_forward(d_vel, d_gp0, d_gp1, d_gp2, d_gvx, d_gvz, d_convvx, d_convvz, d_convpx, d_convpz, d_bx1, d_bz1, d_bx2, d_bz2);
			ptr=d_gp0; d_gp0=d_gp1; d_gp1=d_gp2; d_gp2=ptr;


			cuda_cross_correlate<<<dimg0, dimb0>>>(d_Isg, d_Iss, d_sp0, d_gp0, npml, nnz, nnx);
		}
		cuda_imaging<<<dimg0, dimb0>>>(d_Isg, d_Iss, d_I1, d_I2, npml, nnz, nnx);

		cudaEventRecord(stop);
  		cudaEventSynchronize(stop);
  		cudaEventElapsedTime(&mstimer, start, stop);
    		sf_warning("%d shot finished: %g (s)",is+1, mstimer*1.e-3);
    	}
	cuda_laplace_filter<<<dimg0,dimb0>>>(d_I1,d_sp0,_dz,_dx, npml, nnz, nnx);
	cudaMemcpy(p, d_sp0, N*sizeof(float), cudaMemcpyDeviceToHost);
	window(v0, p, npml, nnz, nnx, nz1, nx1);
	sf_floatwrite(v0,nz1*nx1,imag1); 
	cuda_laplace_filter<<<dimg0,dimb0>>>(d_I2,d_gp0,_dz,_dx, npml, nnz, nnx);
	cudaMemcpy(p, d_gp0, N*sizeof(float), cudaMemcpyDeviceToHost);
	window(v0, p, npml, nnz, nnx, nz1, nx1);
	sf_floatwrite(v0,nz1*nx1,imag2); 
/*
	cudaMemcpy(p, d_I1, N*sizeof(float), cudaMemcpyDeviceToHost);
	window(v0, p, npml, nnz, nnx, nz1, nx1);
	sf_floatwrite(v0,nz1*nx1,imag1); 

	cudaMemcpy(p, d_I2, N*sizeof(float), cudaMemcpyDeviceToHost);
	window(v0, p, npml, nnz, nnx, nz1, nx1);
	sf_floatwrite(v0,nz1*nx1,imag2); 
*/
	cudaEventDestroy(start);
	cudaEventDestroy(stop);

    	free(seis);
    	free(v0);
    	free(vel);
    	free(p);
	device_free();

	return 0;
}
