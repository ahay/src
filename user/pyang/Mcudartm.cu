/* CUDA RTM
 Note: it is not recommended to compile this program using Madgascar.
 The main reason is that Madgascar does a bad job to compile a fast 
	CUDA implementation. The speed will be degraded a lot.
*/
/*
  Copyright (C) 2013  Xi'an Jiaotong University (Pengliang Yang)

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


extern "C" {
#include <rsf.h>
}

#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef true
#define true    (1)
#endif
#ifndef false
#define false   (0)
#endif
#ifndef MAX
#define	MAX(x,y) ((x) > (y) ? (x) : (y))
#endif
#ifndef MIN
#define	MIN(x,y) ((x) < (y) ? (x) : (y))
#endif
#ifndef PI
#define PI (3.141592653589793)
#endif


#include <cuda.h>
#include <cuda_runtime_api.h>

#define Block_Size1 16	// 1st dim of seismic data
#define Block_Size2 16	// 2nd dim of seismic data
const int npml=32;
const int NJ=6;		// finite difference order: NJ=2*k
const int nbell=2;	// radius of Gaussian bell 
const bool frsf=true;	// if true, free surface boundary condition

#include "rtm_kernels.cu"

const bool csdgather=true; // common shot gather (CSD) or not 
char *model_file="true_model.sgy";
const float 	dx=5.0;
const float 	dz=5.0;
const float 	fm=25.0;
const float 	dt=0.001;
const int 	nt=1900;
const int 	ns=20;
const int 	ng=60;

int jsx=13;
int jsz=0;
int jgx=1;
int jgz=0;
int sxbeg=0;//x-begin point of source, index starting from 0
int szbeg=1;//z-begin point of source, index starting from 0
int gxbeg=0;//x-begin point of geophone, index starting from 0
int gzbeg=2;//z-begin point of geophone, index starting from 0


const float 	vmute=1550;
const int 	ntd=int(1.0/(fm*dt));
const float	_dx=1.0/dx;
const float	_dz=1.0/dz;
const int 	nt_h=0.65*nt; // 65% points allocated on host using zero-copy pinned memory; the rest on device

static int 	nz1, nx1, nz, nx, nnz, nnx, N;
static dim3 	dimbbell, dimg0, dimb0;
static dim3 	dimglr1, dimblr1, dimglr2, dimblr2;//lr=left and right
static dim3 	dimgtb1, dimbtb1, dimgtb2, dimbtb2;//tb=top and bottom

// variables on host
float 	*seis, *v0, *vel, *p;
// variables on device
int 	*d_Sxz, *d_Gxz;				// set source and geophone position
float 	*d_bell,*d_wlt, *d_seis,  *d_vel;	// bell, wavelet, seismograms, velocity (vel)
float 	*d_usp,*d_sp0, *d_sp1, *d_svx, *d_svz;		// p, vx, vz for sources
float 	*d_ugp,*d_gp0, *d_gp1, *d_gvx, *d_gvz;			// p, vx, vz for geophones
float 	*d_bx1, *d_bx2, *d_bz1, *d_bz2;		// PML ABC coefficients for p and v (vx, vz)
float 	*d_convpx, *d_convpz, *d_convvx, *d_convvz;// auxiliary variables to decay p and v in PML zone
float 	*d_Iss, *d_Isg, *d_I;			// image result
float 	*h_boundary, *d_boundary;		// boundary on host and device

float	*ptr=NULL;

void matrix_transpose(float *matrix, int nx, int nz)
{
	float *tmp=(float*)malloc(nx*nz*sizeof(float));
	if (tmp==NULL) {printf("out of memory!"); exit(1);}
	for(int iz=0; iz<nz; iz++)
		for(int ix=0; ix<nx; ix++)
			tmp[iz+nz*ix]=matrix[ix+nx*iz];
	memcpy(matrix, tmp, nx*nz*sizeof(float));
	free(tmp);
}

/*
 2nd-order stability condition:	min(dx, dz)>sqrt(2)*dt*max(v)
 numerical dispersion condition: max(dx, dz)<min(v)/(5*fmax)
*/
void check_gird_sanity()
{
	float C;
	if(NJ==2) 		C=1;
	else if (NJ==4)	 	C=0.857;
	else if (NJ==6)		C=0.8;
	else if (NJ==8) 	C=0.777;
	else if (NJ==10)	C=0.759;

	float maxvel=vel[0], minvel=vel[0];
	for(int i=0; i<N; i++)
	{
		if(vel[i]>maxvel) maxvel=vel[i];
		if(vel[i]<minvel) minvel=vel[i];
	}
	float tmp=dt*maxvel*sqrt(1.0/(dx*dx)+1.0/(dz*dz));

	if (tmp>=C) printf("Stability condition not satisfied!\n");
	if (fm>=minvel/(5.0*max(dx,dz))) printf("Dispersion relation not satisfied!\n");
}


void device_alloc()
{
	cudaMalloc(&d_bell,	(2*nbell+1)*(2*nbell+1)*sizeof(float));
    	cudaMalloc(&d_Sxz,	ns*sizeof(int));
    	cudaMalloc(&d_Gxz, 	ng*sizeof(int));
	cudaMalloc(&d_wlt,	(nt+ntd)*sizeof(float));
    	cudaMalloc(&d_seis, 	ng*nt*sizeof(float));
    	cudaMalloc(&d_vel, 	N*sizeof(float));
    	cudaMalloc(&d_usp, 	N*sizeof(float));
    	cudaMalloc(&d_sp0, 	N*sizeof(float));
    	cudaMalloc(&d_sp1, 	N*sizeof(float));
    	cudaMalloc(&d_svx, 	N*sizeof(float));
    	cudaMalloc(&d_svz, 	N*sizeof(float));
    	cudaMalloc(&d_ugp, 	N*sizeof(float));
    	cudaMalloc(&d_gp0, 	N*sizeof(float));
    	cudaMalloc(&d_gp1, 	N*sizeof(float));
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
    	cudaMalloc(&d_I, 	N*sizeof(float));

	cudaHostAlloc(&h_boundary, nt_h*2*(NJ-1)*(nx+nz)*sizeof(float), cudaHostAllocMapped);	
	cudaMalloc(&d_boundary, (nt-nt_h)*2*(NJ-1)*(nx+nz)*sizeof(float));

	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err) 
	printf("Cuda error: Failed to allocate required memory! %s\n", cudaGetErrorString(err));
}


void device_free()
{
	cudaFree(d_bell);
    	cudaFree(d_Sxz);
    	cudaFree(d_Gxz);
	cudaFree(d_wlt);
    	cudaFree(d_seis);
    	cudaFree(d_vel);
    	cudaFree(d_usp);
    	cudaFree(d_sp0);
    	cudaFree(d_sp1);
    	cudaFree(d_svx);
    	cudaFree(d_svz);
    	cudaFree(d_ugp);
    	cudaFree(d_gp0);
    	cudaFree(d_gp1);
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
	cudaFree(d_I);

	cudaFreeHost(h_boundary);
    	cudaFree(d_boundary);

	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err) 
	printf("Cuda error: Failed to free the allocated memory! %s\n", cudaGetErrorString(err));
}

//up-->laplace(p)
void wavefield_ini(float *d_up, float *d_p0, float *d_p1, float *d_vx, float *d_vz, float *d_convpx, float *d_convpz, float *d_convvx, float *d_convvz)
{
	cudaMemset(d_up, 	0,  	N*sizeof(float));
	cudaMemset(d_p0, 	0,	N*sizeof(float));
	cudaMemset(d_p1, 	0,	N*sizeof(float));
	cudaMemset(d_vx, 	0,	N*sizeof(float));
	cudaMemset(d_vz, 	0,	N*sizeof(float));
	cudaMemset(d_convpx, 	0,	2*npml*nnz*sizeof(float));
	cudaMemset(d_convpz, 	0,	2*npml*nnx*sizeof(float));
	cudaMemset(d_convvx, 	0,	2*npml*nnz*sizeof(float));
	cudaMemset(d_convvz, 	0,	2*npml*nnx*sizeof(float));

	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err) 
	printf("Cuda error: Failed to initialize the wavefield variables! %s\n", cudaGetErrorString(err));
}

void forward_laplacian(float *d_up, float *d_p1, float *d_vx, float *d_vz, float *d_convvx, float *d_convvz, float *d_convpx, float *d_convpz, float *d_bx1, float *d_bz1, float *d_bx2, float *d_bz2)
{
	// p0: p{it-1}; p1: p{it};
	if (NJ==2)
	{
		cuda_forward_v_2<<<dimg0, dimb0>>>(d_p1, d_vx, d_vz, _dx, _dz, npml, nnz, nnx);
		cuda_PML_vz_2<<<dimgtb1, dimbtb1>>>(d_p1, d_convpz, d_bz2, d_vz, _dz, npml, nnz, nnx);
		cuda_PML_vx_2<<<dimglr1, dimblr1>>>(d_p1, d_convpx, d_bx2, d_vx, _dx, npml, nnz, nnx);
		cuda_forward_up_2<<<dimg0, dimb0>>>(d_up, d_vx, d_vz, _dx, _dz, npml, nnz, nnx);
		cuda_PML_upz_2<<<dimgtb1, dimbtb1>>>(d_up, d_convvz, d_bz1, d_vz, _dz, npml, nnz, nnx);
		cuda_PML_upx_2<<<dimglr1, dimblr1>>>(d_up, d_convvx, d_bx1, d_vx, _dx, npml, nnz, nnx);
	}
	else if (NJ==4)
	{
		cuda_forward_v_4<<<dimg0, dimb0>>>(d_p1, d_vx, d_vz, _dx, _dz, npml, nnz, nnx);
		cuda_PML_vz_4<<<dimgtb1, dimbtb1>>>(d_p1, d_convpz, d_bz2, d_vz, _dz, npml, nnz, nnx);
		cuda_PML_vx_4<<<dimglr1, dimblr1>>>(d_p1, d_convpx, d_bx2, d_vx, _dx, npml, nnz, nnx);
		cuda_forward_up_4<<<dimg0, dimb0>>>(d_up, d_vx, d_vz, _dx, _dz, npml, nnz, nnx);
		cuda_PML_upz_4<<<dimgtb1, dimbtb1>>>(d_up, d_convvz, d_bz1, d_vz, _dz, npml, nnz, nnx);
		cuda_PML_upx_4<<<dimglr1, dimblr1>>>(d_up, d_convvx, d_bx1, d_vx, _dx, npml, nnz, nnx);
	}
	else if (NJ==6)
	{
		cuda_forward_v_6<<<dimg0, dimb0>>>(d_p1, d_vx, d_vz, _dx, _dz, npml, nnz, nnx);
		cuda_PML_vz_6<<<dimgtb1, dimbtb1>>>(d_p1, d_convpz, d_bz2, d_vz, _dz, npml, nnz, nnx);
		cuda_PML_vx_6<<<dimglr1, dimblr1>>>(d_p1, d_convpx, d_bx2, d_vx, _dx, npml, nnz, nnx);
		cuda_forward_up_6<<<dimg0, dimb0>>>(d_up, d_vx, d_vz, _dx, _dz, npml, nnz, nnx);
		cuda_PML_upz_6<<<dimgtb1, dimbtb1>>>(d_up, d_convvz, d_bz1, d_vz, _dz, npml, nnz, nnx);
		cuda_PML_upx_6<<<dimglr1, dimblr1>>>(d_up, d_convvx, d_bx1, d_vx, _dx, npml, nnz, nnx);
	}
	else if (NJ==8)
	{
		cuda_forward_v_8<<<dimg0, dimb0>>>(d_p1, d_vx, d_vz, _dx, _dz, npml, nnz, nnx);
		cuda_PML_vz_8<<<dimgtb1, dimbtb1>>>(d_p1, d_convpz, d_bz2, d_vz, _dz, npml, nnz, nnx);
		cuda_PML_vx_8<<<dimglr1, dimblr1>>>(d_p1, d_convpx, d_bx2, d_vx, _dx, npml, nnz, nnx);
		cuda_forward_up_8<<<dimg0, dimb0>>>(d_up, d_vx, d_vz, _dx, _dz, npml, nnz, nnx);
		cuda_PML_upz_8<<<dimgtb1, dimbtb1>>>(d_up, d_convvz, d_bz1, d_vz, _dz, npml, nnz, nnx);
		cuda_PML_upx_8<<<dimglr1, dimblr1>>>(d_up, d_convvx, d_bx1, d_vx, _dx, npml, nnz, nnx);
	}
	else if (NJ==10)
	{
		cuda_forward_v_10<<<dimg0, dimb0>>>(d_p1, d_vx, d_vz, _dx, _dz, npml, nnz, nnx);
		cuda_PML_vz_10<<<dimgtb1, dimbtb1>>>(d_p1, d_convpz, d_bz2, d_vz, _dz, npml, nnz, nnx);
		cuda_PML_vx_10<<<dimglr1, dimblr1>>>(d_p1, d_convpx, d_bx2, d_vx, _dx, npml, nnz, nnx);
		cuda_forward_up_10<<<dimg0, dimb0>>>(d_up, d_vx, d_vz, _dx, _dz, npml, nnz, nnx);
		cuda_PML_upz_10<<<dimgtb1, dimbtb1>>>(d_up, d_convvz, d_bz1, d_vz, _dz, npml, nnz, nnx);
		cuda_PML_upx_10<<<dimglr1, dimblr1>>>(d_up, d_convvx, d_bx1, d_vx, _dx, npml, nnz, nnx);
	}
}

void backward_laplacian(float *d_up, float *d_p1, float *d_vx, float *d_vz)
{
	// p0: p{it-1}; p1: p{it};
	if (NJ==2)
	{
		cuda_forward_v_2<<<dimg0, dimb0>>>(d_p1, d_vx, d_vz, _dx, _dz, npml, nnz, nnx);
		cuda_forward_up_2<<<dimg0, dimb0>>>(d_up, d_vx, d_vz, _dx, _dz, npml, nnz, nnx);
	}
	else if (NJ==4)
	{
		cuda_forward_v_4<<<dimg0, dimb0>>>(d_p1, d_vx, d_vz, _dx, _dz, npml, nnz, nnx);
		cuda_forward_up_4<<<dimg0, dimb0>>>(d_up, d_vx, d_vz, _dx, _dz, npml, nnz, nnx);
	}
	else if (NJ==6)
	{
		cuda_forward_v_6<<<dimg0, dimb0>>>(d_p1, d_vx, d_vz, _dx, _dz, npml, nnz, nnx);
		cuda_forward_up_6<<<dimg0, dimb0>>>(d_up, d_vx, d_vz, _dx, _dz, npml, nnz, nnx);
	}
	else if (NJ==8)
	{
		cuda_forward_v_8<<<dimg0, dimb0>>>(d_p1, d_vx, d_vz, _dx, _dz, npml, nnz, nnx);
		cuda_forward_up_8<<<dimg0, dimb0>>>(d_up, d_vx, d_vz, _dx, _dz, npml, nnz, nnx);
	}
	else if (NJ==10)
	{
		cuda_forward_v_10<<<dimg0, dimb0>>>(d_p1, d_vx, d_vz, _dx, _dz, npml, nnz, nnx);
		cuda_forward_up_10<<<dimg0, dimb0>>>(d_up, d_vx, d_vz, _dx, _dz, npml, nnz, nnx);
	}
}


int main(int argc, char *argv[])
{
	nz1=320;
	nx1=320;
	nx=(int)((nx1+Block_Size1-1)/Block_Size1)*Block_Size1;
	nz=(int)((nz1+Block_Size2-1)/Block_Size2)*Block_Size2;
    	nnz = 2*npml+nz;
    	nnx = 2*npml+nx;
	N=nnz*nnx;
	dimbbell=dim3(2*nbell+1,2*nbell+1);
    	dimb0=dim3(Block_Size1, Block_Size2);  	dimg0=dim3(nnz/Block_Size1, nnx/Block_Size2);
	dimblr1=dim3(Block_Size1, 32);		dimglr1=dim3(nnz/Block_Size1, 2);
	dimbtb1=dim3(32, Block_Size2); 		dimgtb1=dim3(2, nnx/Block_Size2);
	dimblr2=dim3(nz/Block_Size1,(NJ+15)/16); dimglr2=dim3(Block_Size1, 16);
	dimbtb2=dim3(16, Block_Size2);		dimgtb2=dim3((NJ+15)/16, nx/Block_Size2);

	
    	seis=(float*)malloc(ng*nt*sizeof(float));
	if (seis==NULL) { printf("out of memory!"); exit(1);}
    	v0=(float*)malloc(nx*nz*sizeof(float));
	if (v0==NULL) 	{ printf("out of memory!"); exit(1);}
    	vel=(float*)malloc(N*sizeof(float));
	if (vel==NULL) 	{ printf("out of memory!"); exit(1);}
    	p=(float*)malloc(N*sizeof(float));
	if (p==NULL) 	{ printf("out of memory!"); exit(1);}
    	memset(v0, 0, nz1*nx1*sizeof(float));
    	memset(vel, 0, N*sizeof(float));
    	memset(p, 0, N*sizeof(float));

	printf("NJ=%d \t\n",	NJ);
    	printf("npml=%d \t\n",	npml);
    	printf("nx1=%d   \t\n",	nx1);
    	printf("nz1=%d   \t\n",	nz1);
    	printf("nx=%d   \t\n",	nx);
    	printf("nz=%d   \t\n",	nz);
    	printf("nnx=%d  \t\n",	nnx);
    	printf("nnz=%d  \t\n",	nnz);
    	printf("dx=%g\t\t(m)\n", dx);
    	printf("dz=%g\t\t(m)\n", dz);
    	printf("dt=%g  \t(s)\n", dt);
    	printf("fm=%g\t\t(Hz)\n", fm);
    	printf("nt=%d\n",	nt);
   	printf("ns=%d\n",	ns);
    	printf("ng=%d\n",	ng);
	check_gird_sanity();


    	cudaSetDevice (0);
	device_alloc();

	cuda_ini_bell<<<dim3(1,1),dim3(2*nbell+1,2*nbell+1)>>>(d_bell);
	cuda_ricker_wavelet<<<(nt+ntd+511)/512,512>>>(d_wlt, fm, dt, nt+ntd);

	if (!(sxbeg>=0 && szbeg>=0 && sxbeg+(ns-1)*jsx<nx && szbeg+(ns-1)*jsz<nz))	
	{ printf("sources exceeds the computing zone!\n"); exit(1);}
	cuda_set_sg<<<(ns+255)/256, 256>>>(d_Sxz, sxbeg, szbeg, jsx, jsz, ns, npml, nnz);	


	int distx=sxbeg-gxbeg;
	int distz=szbeg-gzbeg;
	if (csdgather)
	{
		//distance between source and geophone at the beginning
		if (!(gxbeg>=0 && gzbeg>=0 && gxbeg+(ng-1)*jgx<nx && gzbeg+(ng-1)*jgz<nz &&
		(sxbeg+(ns-1)*jsx)+(ng-1)*jgx-distx <nx  && (szbeg+(ns-1)*jsz)+(ng-1)*jgz-distz <nz))	
		{ printf("geophones exceeds the computing zone!\n"); exit(1);}
		cuda_set_sg<<<(ng+255)/256, 256>>>(d_Gxz, gxbeg, gzbeg, jgx, jgz, ng, npml, nnz);
	}
	else{
		if (!(gxbeg>=0 && gzbeg>=0 && gxbeg+(ng-1)*jgx<nx && gzbeg+(ng-1)*jgz<nz))	
		{ printf("geophones exceeds the computing zone!\n"); exit(1);}
		cuda_set_sg<<<(ng+255)/256, 256>>>(d_Gxz, gxbeg, gzbeg, jgx, jgz, ng, npml, nnz);
	}

	cuda_set2<<<dimg0, dimb0>>>(d_vel, 1800.0, 2000.0,nnz, nnx);
    	//cudaMemcpy(d_vel, vel, N*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemset(d_Iss, 0,  	N*sizeof(float));
	cudaMemset(d_Isg, 0,  	N*sizeof(float));
	cudaMemset(d_I, 0,  	N*sizeof(float));
    	

	float mstimer = 0;// timer unit: millionseconds
	cudaEvent_t start, stop;
  	cudaEventCreate(&start);	
	cudaEventCreate(&stop);

	FILE *fp=fopen("wav.dat","wb");
	if (fp==NULL) { printf("cannot open the file"); exit(1);}
	cuda_ini_abcz<<<dimgtb1, dimbtb1>>>(d_vel, d_bz1, d_bz2, dx, dz, dt, npml, nnz, nnx);
	cuda_ini_abcx<<<dimglr1, dimblr1>>>(d_vel, d_bx1, d_bx2, dx, dz, dt, npml, nnz, nnx);
    	for(int is=0; is<ns; is++)
    	{
		cudaEventRecord(start);

		cudaMemset(d_Iss, 	0,	N*sizeof(float));
	    	cudaMemset(d_Isg, 	0,	N*sizeof(float));
		cudaMemset(d_seis, 	0,	nt*ng*sizeof(float));
		cudaMemset(h_boundary, 	0,	nt_h*2*(NJ-1)*(nx+nz)*sizeof(float));
	    	cudaMemset(d_boundary, 	0,	(nt-nt_h)*2*(NJ-1)*(nx+nz)*sizeof(float));
		wavefield_ini(d_usp, d_sp0, d_sp1, d_svx, d_svz, d_convpx, d_convpz, d_convvx, d_convvz);
		if (csdgather)
		{
			gxbeg=sxbeg+is*jsx-distx;
			cuda_set_sg<<<(ng+255)/256, 256>>>(d_Gxz, gxbeg, gzbeg, jgx, jgz, ng, npml, nnz);
		}
		for(int it=0; it<nt+ntd; it++)
		{
			int kt=it-ntd;
			forward_laplacian(d_usp, d_sp1, d_svx, d_svz, d_convvx, d_convvz, d_convpx, d_convpz, d_bx1, d_bz1, d_bx2, d_bz2);
			cuda_add_bellwlt<<<dim3(1,1), dimbbell>>>(d_usp, d_bell, &d_wlt[it], &d_Sxz[is], 1, npml, nnz, nnx, true);
			cuda_step_forward<<<dimg0,dimb0>>>(d_vel, d_usp, d_sp0, d_sp1, dt, npml, nnz, nnx);
			ptr=d_sp0; d_sp0=d_sp1; d_sp1=ptr;

			if (kt>=0)
			{
				cuda_record<<<(ng+255)/256, 256>>>(d_sp1, &d_seis[kt*ng], d_Gxz, ng);
/*
				if (kt%50==0)
				{
					cudaMemcpy(p, d_sp1, N*sizeof(float), cudaMemcpyDeviceToHost);
					fwrite(p, sizeof(float), N, fp);
				}
*/
				cuda_mute<<<(ng+511)/512, 512>>>(&d_seis[kt*ng], gzbeg, szbeg, gxbeg, sxbeg+is*jsx, jgx, kt, ntd, vmute, dt, dz, dx, ng);

				if(kt<nt_h) cudaHostGetDevicePointer(&ptr, &h_boundary[kt*2*(NJ-1)*(nx+nz)], 0);
				else  ptr=&d_boundary[(kt-nt_h)*2*(NJ-1)*(nx+nz)];
				cuda_rw_innertb<<<dimgtb2, dimbtb2>>>(ptr, 	d_sp0, npml, nnz, nnx, NJ, false);
				cuda_rw_innerlr<<<dimglr2, dimblr2>>>(&ptr[2*(NJ-1)*nx], d_sp0, npml, nnz, nnx, NJ, false);
			}		    
		}
/*
		cudaMemcpy(seis, d_seis, nt*ng*sizeof(float), cudaMemcpyDeviceToHost);
		matrix_transpose(seis, ng, nt);// before: nx=ng; nz=nt; after: nx=nt; nz=ng;
		fwrite(seis, sizeof(float), nt*ng, fp);
*/

		ptr=d_sp0; d_sp0=d_sp1; d_sp1=ptr;
		// revese time loop
		wavefield_ini(d_ugp,d_gp0, d_gp1, d_gvx, d_gvz, d_convpx, d_convpz, d_convvx, d_convvz);
		for(int it=nt-1+ntd; it>-1+ntd; it--)
		{
			int kt=it-ntd;
			// backward time step receiver wavefield
			forward_laplacian(d_ugp, d_gp1, d_gvx, d_gvz, d_convvx, d_convvz, d_convpx, d_convpz, d_bx1, d_bz1, d_bx2, d_bz2);
			// add receiver term
			cuda_add_wavelet<<<(ng+255)/256,256>>>(d_ugp, &d_seis[kt*ng], d_Gxz, ng, true);
			cuda_step_forward<<<dimg0,dimb0>>>(d_vel, d_ugp, d_gp0, d_gp1, dt,npml, nnz, nnx);
			ptr=d_gp0; d_gp0=d_gp1; d_gp1=ptr;

			// read saved boundary
			if(kt<nt_h) cudaHostGetDevicePointer(&ptr, &h_boundary[kt*2*(NJ-1)*(nx+nz)], 0);
			else  ptr=&d_boundary[(kt-nt_h)*2*(NJ-1)*(nx+nz)];
			cuda_rw_innertb<<<dimgtb2, dimbtb2>>>(ptr, 		d_sp1, npml, nnz, nnx, NJ, true);
			cuda_rw_innerlr<<<dimglr2, dimblr2>>>(&ptr[2*(NJ-1)*nx], d_sp1, npml, nnz, nnx, NJ, true);
			backward_laplacian(d_usp, d_sp1, d_svx, d_svz);
			// subtract the wavelet
			cuda_add_bellwlt<<<dim3(1,1), dimbbell>>>(d_usp, d_bell, &d_wlt[it], &d_Sxz[is], 1, npml, nnz, nnx, true);
			// backward time step source wavefield
			cuda_step_forward<<<dimg0,dimb0>>>(d_vel, d_usp, d_sp0, d_sp1, dt, npml, nnz, nnx);
			ptr=d_sp0; d_sp0=d_sp1; d_sp1=ptr;

			// apply imaging condition
			cuda_cross_correlate<<<dimg0, dimb0>>>(d_Isg, d_Iss, d_sp0, d_gp0, npml, nnz, nnx);
		}
		cuda_normalize_image<<<dimg0, dimb0>>>(d_Isg, d_Iss, d_I, npml, nnz, nnx);

		cudaEventRecord(stop);
  		cudaEventSynchronize(stop);
  		cudaEventElapsedTime(&mstimer, start, stop);
    		printf("shot %d  finished: %f (s)\n",is+1, mstimer*1e-3);
    	}

	cudaMemcpy(p, d_I, N*sizeof(float), cudaMemcpyDeviceToHost);
	fwrite(p, sizeof(float), N, fp);
	fclose(fp);

	cudaEventDestroy(start);
	cudaEventDestroy(stop);

    	free(seis);
    	free(v0);
    	free(vel);
    	free(p);

	device_free();

	
	exit(0);
}
