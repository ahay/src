/* Generate shots for FWI using Enquist absorbing boundary condition

Note: 	You can try other complex boundary condition but we do not
	recommend to do so. The main reason is that FWI is to recover
	the low-frequency information of the earth model. Low-freq 
	means that exact absorbing is not necessarilly needed. The 
	result will be improved with the optimization procedure. 
	Furthermore, complex boundary condition (such as sponge ABC or
	PML) implies more computational cost, which will degrade the
	efficiency of FWI. 
*/
/*
  Copyright (C) 2013  Xi'an Jiaotong University (Pengliang Yang)
    Email: ypl.2100@gmail.com	
    Acknowledgement: This code is written with the help of Baoli Wang.

References:
    [1] Clayton, Robert, and Björn Engquist. "Absorbing boundary 
	conditions for acoustic and elastic wave equations." Bulletin 
	of the Seismological Society of America 67.6 (1977): 1529-1540.
    [2] Tarantola, Albert. "Inversion of seismic reflection data in the 
	acoustic approximation." Geophysics 49.8 (1984): 1259-1266.
    [3] Pica, A., J. P. Diet, and A. Tarantola. "Nonlinear inversion 
	of seismic reflection data in a laterally invariant medium." 
	Geophysics 55.3 (1990): 284-292.
    [4] Hager, William W., and Hongchao Zhang. "A survey of nonlinear
	conjugate gradient methods." Pacific journal of Optimization 
	2.1 (2006): 35-58.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <cuda_runtime.h>

extern "C" {
#include <rsf.h>
}

#ifndef true
#define true    (1)
#endif
#ifndef false
#define false   (0)
#endif
#ifndef EPS
#define EPS	1.0e-15f
#endif

#define PI 	3.141592653589793f
#define Block_Size1 16	/* 1st dim block size */
#define Block_Size2 16  /* 2nd dim block size */
#define Block_Size  512	/* vector computation blocklength */
#define nbell	2	/* radius of Gaussian bell: diameter=2*nbell+1 */

#include "cuda_fwi_kernels.cu"

static bool csdgather;
static int nz,nx,nz1,nx1,nt,ns,ng;
static float dx, dz, fm, dt;


/* variables on host */
float 	*v0, *dobs, *vv;
/* variables on device */
int 	*d_sxz, *d_gxz;			
float 	*d_wlt, *d_vv, *d_sp0, *d_sp1, *d_dobs;

void matrix_transpose(float *matrix, float *trans, int n1, int n2)
/*< matrix transpose: matrix tansposed to be trans >*/
{
	int i1, i2;

	for(i2=0; i2<n2; i2++)
	for(i1=0; i1<n1; i1++)
	    trans[i2+n2*i1]=matrix[i1+n1*i2];
}

void expand(float*vv, float *v0, int nz, int nx, int nz1, int nx1)
/*< round up the model size to be multiples of block size >*/
{
	int i1,i2,i11,i22;

	for(i2=0; i2<nx; i2++)
	for(i1=0; i1<nz; i1++)
	{
		i11=(i1<nz1)?i1:(nz1-1);
		i22=(i2<nx1)?i2:(nx1-1);
		vv[i1+i2*nz]=v0[i11+nz1*i22];
	}	
}


void window(float *v0,float *vv, int nz, int nx, int nz1, int nx1)
/*< window the portion to be the same size as initial model >*/
{
	int i1, i2;

	for(i2=0; i2<nx1; i2++)
	for(i1=0; i1<nz1; i1++)
		  v0[i1+i2*nz1]=vv[i1+nz*i2];
}

void device_alloc()
/*< allocate memories for variables on device >*/
{
	cudaMalloc(&d_vv, nz*nx*sizeof(float));
	cudaMalloc(&d_sp0, nz*nx*sizeof(float));
	cudaMalloc(&d_sp1, nz*nx*sizeof(float));
	cudaMalloc(&d_wlt, nt*sizeof(float));
	cudaMalloc(&d_sxz, nt*sizeof(float));
	cudaMalloc(&d_gxz, ng*sizeof(float));
	cudaMalloc(&d_dobs, ng*nt*sizeof(float));

    	cudaError_t err = cudaGetLastError ();
    	if (cudaSuccess != err) 
	printf("Cuda error: Failed to allocate required memory!: %s\n", cudaGetErrorString(err));
}


void device_free()
/*< free the variables on device >*/
{
	cudaFree(d_vv);
	cudaFree(d_sp0);
	cudaFree(d_sp1);
	cudaFree(d_wlt);
	cudaFree(d_sxz);
	cudaFree(d_gxz);
	cudaFree(d_dobs);

    	cudaError_t err = cudaGetLastError ();
    	if (cudaSuccess != err)
	printf("Cuda error: Failed to free the allocated memory!: %s\n", cudaGetErrorString(err));
}


int main(int argc, char *argv[])
{
	int is, it, distx, distz,sxbeg,szbeg,gxbeg,gzbeg,jsx,jsz,jgx,jgz;
	float dtx,dtz,mstimer,amp;
	float *trans, *ptr=NULL;
	sf_file vinit, shots;

    	/* initialize Madagascar */
    	sf_init(argc,argv);

    	/*< set up I/O files >*/
    	vinit=sf_input ("in");   /* initial velocity model, unit=m/s */
    	shots=sf_output("out");  /* output image with correlation imaging condition */ 

    	/* get parameters for forward modeling */
    	if (!sf_histint(vinit,"n1",&nz1)) sf_error("no n1");
    	if (!sf_histint(vinit,"n2",&nx1)) sf_error("no n2");
    	if (!sf_histfloat(vinit,"d1",&dz)) sf_error("no d1");
   	if (!sf_histfloat(vinit,"d2",&dx)) sf_error("no d2");

	if (!sf_getfloat("amp",&amp)) amp=1000;
	/* maximum amplitude of ricker */
    	if (!sf_getfloat("fm",&fm)) fm=10;	
	/* dominant freq of ricker */
    	if (!sf_getfloat("dt",&dt)) sf_error("no dt");	
	/* time interval */
    	if (!sf_getint("nt",&nt))   sf_error("no nt");	
	/* total modeling time steps */
    	if (!sf_getint("ns",&ns))   sf_error("no ns");	
	/* total shots */
    	if (!sf_getint("ng",&ng))   sf_error("no ng");	
	/* total receivers in each shot */	
    	if (!sf_getint("jsx",&jsx))   sf_error("no jsx");
	/* source x-axis  jump interval  */
    	if (!sf_getint("jsz",&jsz))   jsz=0;
	/* source z-axis jump interval  */
    	if (!sf_getint("jgx",&jgx))   jgx=1;
	/* receiver x-axis jump interval */
    	if (!sf_getint("jgz",&jgz))   jgz=0;
	/* receiver z-axis jump interval */
    	if (!sf_getint("sxbeg",&sxbeg))   sf_error("no sxbeg");
	/* x-begining index of sources, starting from 0 */
    	if (!sf_getint("szbeg",&szbeg))   sf_error("no szbeg");
	/* z-begining index of sources, starting from 0 */
    	if (!sf_getint("gxbeg",&gxbeg))   sf_error("no gxbeg");
	/* x-begining index of receivers, starting from 0 */
    	if (!sf_getint("gzbeg",&gzbeg))   sf_error("no gzbeg");
	/* z-begining index of receivers, starting from 0 */
	if (!sf_getbool("csdgather",&csdgather)) csdgather=false;
	/* default, common shot-gather; if n, record at every point*/

	sf_putint(shots,"n1",nt);	
	sf_putint(shots,"n2",ng);
	sf_putint(shots,"n3",ns);
	sf_putfloat(shots,"d1",dt);
	sf_putfloat(shots,"d2",jgx*dx);
	sf_putfloat(shots,"o1",0);
	sf_putstring(shots,"label1","Time");
	sf_putstring(shots,"label2","Lateral");
	sf_putstring(shots,"label3","Shot");
	sf_putstring(shots,"unit1","sec");
	sf_putstring(shots,"unit2","m");
	sf_putfloat(shots,"amp",amp);
	sf_putfloat(shots,"fm",fm);
	sf_putint(shots,"ng",ng);
	sf_putint(shots,"szbeg",szbeg);
	sf_putint(shots,"sxbeg",sxbeg);
	sf_putint(shots,"gzbeg",gzbeg);
	sf_putint(shots,"gxbeg",gxbeg);
	sf_putint(shots,"jsx",jsx);
	sf_putint(shots,"jsz",jsz);
	sf_putint(shots,"jgx",jgx);
	sf_putint(shots,"jgz",jgz);
	sf_putint(shots,"csdgather",csdgather?1:0);

	dtx=dt/dx; 
	dtz=dt/dz; 
	/* round the size up to multiples of Block size */
	nx=(int)((nx1+Block_Size1-1)/Block_Size1)*Block_Size1;
	nz=(int)((nz1+Block_Size2-1)/Block_Size2)*Block_Size2;

	v0=(float*)malloc(nz1*nx1*sizeof(float));
	vv=(float*)malloc(nz*nx*sizeof(float));
	dobs=(float*)malloc(ng*nt*sizeof(float));
	trans=(float*)malloc(ng*nt*sizeof(float));
	sf_floatread(v0,nz1*nx1,vinit);
	expand(vv, v0, nz, nx, nz1, nx1);
	memset(dobs,0,ng*nt*sizeof(float));
	memset(trans,0,ng*nt*sizeof(float));

    	cudaSetDevice(0);
    	cudaError_t err = cudaGetLastError();
    	if (cudaSuccess != err) 
	printf("Cuda error: Failed to initialize device: %s\n", cudaGetErrorString(err));
	device_alloc(); 

	dim3 dimg=dim3(nz/Block_Size1, nx/Block_Size2),dimb=dim3(Block_Size1, Block_Size2); 

	cudaMemcpy(d_vv, vv, nz*nx*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemset(d_sp0,0,nz*nx*sizeof(float));
	cudaMemset(d_sp1,0,nz*nx*sizeof(float));
	cuda_ricker_wavelet<<<(nt+511)/512,512>>>(d_wlt, amp, fm, dt, nt);
	if (!(sxbeg>=0 && szbeg>=0 && sxbeg+(ns-1)*jsx<nx1 && szbeg+(ns-1)*jsz<nz1))	
	{ printf("sources exceeds the computing zone!\n"); exit(1);}
	cuda_set_sg<<<(ns+511)/512,512>>>(d_sxz, sxbeg, szbeg, jsx, jsz, ns, nz);
	distx=sxbeg-gxbeg;
	distz=szbeg-gzbeg;
	if (csdgather)	{
		if (!(gxbeg>=0 && gzbeg>=0 && gxbeg+(ng-1)*jgx<nx1 && gzbeg+(ng-1)*jgz<nz1 &&
		(sxbeg+(ns-1)*jsx)+(ng-1)*jgx-distx <nx1  && (szbeg+(ns-1)*jsz)+(ng-1)*jgz-distz <nz1))	
		{ printf("geophones exceeds the computing zone!\n"); exit(1);}
	}
	else{
		if (!(gxbeg>=0 && gzbeg>=0 && gxbeg+(ng-1)*jgx<nx1 && gzbeg+(ng-1)*jgz<nz1))	
		{ printf("geophones exceeds the computing zone!\n"); exit(1);}
	}
	cuda_set_sg<<<(ng+511)/512,512>>>(d_gxz, gxbeg, gzbeg, jgx, jgz, ng, nz);

	cudaEvent_t start, stop;
  	cudaEventCreate(&start);	
	cudaEventCreate(&stop);
	for(is=0;is<ns;is++)
	{
		cudaEventRecord(start);
		cudaMemset(d_dobs, 0, ng*nt*sizeof(float));
		if (csdgather)	{
			gxbeg=sxbeg+is*jsx-distx;
			cuda_set_sg<<<(ng+511)/512, 512>>>(d_gxz, gxbeg, gzbeg, jgx, jgz, ng, nz);
		}
		cudaMemset(d_sp0, 0, nz*nx*sizeof(float));
		cudaMemset(d_sp1, 0, nz*nx*sizeof(float));
		for(it=0; it<nt; it++)
		{
			cuda_add_source<<<1,1>>>(d_sp1, &d_wlt[it], &d_sxz[is], 1, true);
			cuda_step_forward<<<dimg,dimb>>>(d_sp0, d_sp1, d_vv, dtz, dtx, nz, nx);
			ptr=d_sp0; d_sp0=d_sp1; d_sp1=ptr;

			cuda_record<<<(ng+511)/512, 512>>>(d_sp0, &d_dobs[it*ng], d_gxz, ng);
		}
		cudaMemcpy(dobs, d_dobs, ng*nt*sizeof(float), cudaMemcpyDeviceToHost);
		matrix_transpose(dobs, trans, ng, nt);
		sf_floatwrite(trans,ng*nt,shots);

		cudaEventRecord(stop);
  		cudaEventSynchronize(stop);
  		cudaEventElapsedTime(&mstimer, start, stop);
    		sf_warning("shot %d finished: %f (s)",is+1, mstimer*1e-3);
	}
	cudaEventDestroy(start);
	cudaEventDestroy(stop);

	free(v0);
	free(vv);
	free(dobs);
	free(trans);
	device_free();

	return 0;
}

