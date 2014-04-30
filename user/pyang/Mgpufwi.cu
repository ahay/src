/* CUDA based FWI using Enquist absorbing boundary condition

Note: 	You can try other complex boundary condition but we do not
	recommend to do so. The main reason is that FWI is to recover
	the low-frequency information of the earth model. Low-freq 
	means that exact absorbing is not necessarilly needed. The 
	result will be improved with the optimization precedure. 
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

#define PI 	3.141592653589793f
#define Block_Size1 16		// 1st dim block size
#define Block_Size2 16		// 2nd dim block size
#define Block_Size  512		// vector computation blocklength
#define nbell	2		// radius of Gaussian bell: diameter=2*nbell+1
const int niter=300;		// total iterations

#include "cuda_fwi_kernels.cu"

/*
const bool csdgather=false; 	// common shot gather (CSD) or not
char *model_file="ini_syn100x120.bin";
char *shots_file="shots.bin";
const int 	nz=100;
const int 	nx=120;
const float 	dx=5.0;
const float 	dz=5.0;
const float 	fm=20.0;
const float 	dt=0.001;
const int 	nt=800;
const int 	ns=12;
const int 	ng=120;

int jsx=10;	//source x-axis jump interval
int jsz=0;	//source z-axis jump interval
int jgx=1;      //geophone x-axis jump interval
int jgz=0;	//geophone z-axis jump interval
int sxbeg=5;	//x-begining index of source, starting from 0
int szbeg=2;	//z-begining index of source, starting from 0
int gxbeg=0;	//x-begining index of geophone, starting from 0
int gzbeg=3;	//z-begining index of geophone, starting from 0
*/

const bool csdgather=false; 	// common shot gather (CSD) or not
char *model_file="sm_marm240x737.bin";
char *shots_file="shots.bin";
const int 	nz=240;
const int 	nx=737;
const float 	dx=12.5;
const float 	dz=12.5;
const float 	fm=15.0;
const float 	dt=0.001;
const int 	nt=3000;
const int 	ns=20;
const int 	ng=737;

int jsx=37;
int jsz=0;
int jgx=1;
int jgz=0;
int sxbeg=15;//x-begin point of source, index starting from 0
int szbeg=2;//z-begin point of source, index starting from 0
int gxbeg=0;//x-begin point of geophone, index starting from 0
int gzbeg=3;//z-begin point of geophone, index starting from 0

// variables on host
float 	*v0, *dobs;
// variables on device
int 	*d_sxz, *d_gxz;			
float 	*d_wlt, *d_vv, *d_sillum, *d_gillum, *d_lap, *d_vtmp, *d_sp0, *d_sp1, *d_gp0, *d_gp1,*d_bndr;
float	*d_dobs, *d_dcal, *d_derr, *d_g0, *d_g1, *d_cg, *d_pars, *d_alpha1, *d_alpha2, *d_bell;
/*
d_pars[0]: obj;
d_pars[1]: beta;
d_pars[2]: epsilon;
d_pars[3]: alpha;
d_alpha1[]: numerator of alpha, length=ng
d_alpha2[]: denominator of alpha, length=ng
*/

void matrix_transpose(float *matrix, int n1, int n2)
/*< matrix transpose >*/
{
	float *tmp=(float*)malloc(n1*n2*sizeof(float));
	if (tmp==NULL) {printf("out of memory!"); exit(1);}
	for(int i2=0; i2<n2; i2++){
		for(int i1=0; i1<n1; i1++){
			tmp[i2+n2*i1]=matrix[i1+n1*i2];
		}
	}
	memcpy(matrix, tmp, n1*n2*sizeof(float));
	free(tmp);
}



void device_alloc()
/*< allocate memories for variables on device >*/
{
	cudaMalloc(&d_vv, nz*nx*sizeof(float));
	cudaMalloc(&d_sp0, nz*nx*sizeof(float));
	cudaMalloc(&d_sp1, nz*nx*sizeof(float));
	cudaMalloc(&d_gp0, nz*nx*sizeof(float));
	cudaMalloc(&d_gp1, nz*nx*sizeof(float));
	cudaMalloc(&d_wlt, nt*sizeof(float));
	cudaMalloc(&d_sxz, nt*sizeof(float));
	cudaMalloc(&d_gxz, ng*sizeof(float));
	cudaMalloc(&d_bndr, nt*2*(nz+nx)*sizeof(float));
	cudaMalloc(&d_dobs, ng*nt*sizeof(float));
	cudaMalloc(&d_dcal, ng*sizeof(float));
	cudaMalloc(&d_derr, ns*ng*nt*sizeof(float));
	cudaMalloc(&d_g0, nz*nx*sizeof(float));
	cudaMalloc(&d_g1, nz*nx*sizeof(float));
	cudaMalloc(&d_cg, nz*nx*sizeof(float));
	cudaMalloc(&d_lap, nz*nx*sizeof(float));
	cudaMalloc(&d_sillum, nz*nx*sizeof(float));
	cudaMalloc(&d_gillum, nz*nx*sizeof(float));
	cudaMalloc(&d_pars, 4*sizeof(float));
	cudaMalloc(&d_alpha1, ng*sizeof(float));
	cudaMalloc(&d_alpha2, ng*sizeof(float));
	cudaMalloc(&d_vtmp, nz*nx*sizeof(float));
	cudaMalloc(&d_bell, (2*nbell+1)*sizeof(float));

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
	cudaFree(d_gp0);
	cudaFree(d_gp1);
	cudaFree(d_wlt);
	cudaFree(d_sxz);
	cudaFree(d_gxz);
	cudaFree(d_bndr);
	cudaFree(d_dobs);
	cudaFree(d_dcal);
	cudaFree(d_derr);
	cudaFree(d_g0);
	cudaFree(d_g1);
	cudaFree(d_cg);
	cudaFree(d_lap);
	cudaFree(d_sillum);
	cudaFree(d_gillum);
	cudaFree(d_pars);
	cudaFree(d_alpha1);
	cudaFree(d_alpha2);
	cudaFree(d_vtmp);
	cudaFree(d_bell);

    	cudaError_t err = cudaGetLastError ();
    	if (cudaSuccess != err)
	printf("Cuda error: Failed to free the allocated memory!: %s\n", cudaGetErrorString(err));
}


void wavefield_init(float *d_p0, float *d_p1, int N)
/*< initialize wavefield >*/
{
	cudaMemset(d_p0, 0, N*sizeof(float));
	cudaMemset(d_p1, 0, N*sizeof(float));

    	cudaError_t err = cudaGetLastError ();
    	if (cudaSuccess != err) 
	printf("Cuda error: Failed to initialize the wavefield variables!: %s\n", cudaGetErrorString(err));
}


void report_par(float *d_par,char *s)
/*< report the parameter with name s[] >*/
{
	float tmp;
	cudaMemcpy(&tmp, d_par, sizeof(float), cudaMemcpyDeviceToHost);
	printf("%s=%f\n",s, tmp);
}


int main(int argc, char *argv[])
{
	v0=(float*)malloc(nz*nx*sizeof(float));
	dobs=(float*)malloc(ng*nt*sizeof(float));
	FILE *fp;
	fp=fopen(model_file,"rb");
	if (fp==NULL) { printf("cannot open file %s:\n",__FILE__); exit(1);}
	fread(v0, sizeof(float), nz*nx, fp);
	fclose(fp);
	memset(dobs,0, ng*nt*sizeof(float));

    	cudaSetDevice(0);
    	cudaError_t err = cudaGetLastError();
    	if (cudaSuccess != err) 
	printf("Cuda error: Failed to initialize device: %s", cudaGetErrorString(err));
	device_alloc(); 


	cudaMemcpy(d_vv, v0, nz*nx*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemset(d_sp0,0,nz*nx*sizeof(float));
	cudaMemset(d_sp1,0,nz*nx*sizeof(float));
	cudaMemset(d_gp0,0,nz*nx*sizeof(float));
	cudaMemset(d_gp1,0,nz*nx*sizeof(float));
	cuda_ricker_wavelet<<<(nt+511)/512,512>>>(d_wlt, fm, dt, nt);
	if (!(sxbeg>=0 && szbeg>=0 && sxbeg+(ns-1)*jsx<nx && szbeg+(ns-1)*jsz<nz))	
	{ printf("sources exceeds the computing zone!\n"); exit(1);}
	cuda_set_sg<<<(ns+511)/512,512>>>(d_sxz, sxbeg, szbeg, jsx, jsz, ns, nz);
	int distx=sxbeg-gxbeg;
	int distz=szbeg-gzbeg;
	if (csdgather)	{
		if (!(gxbeg>=0 && gzbeg>=0 && gxbeg+(ng-1)*jgx<nx && gzbeg+(ng-1)*jgz<nz &&
		(sxbeg+(ns-1)*jsx)+(ng-1)*jgx-distx <nx  && (szbeg+(ns-1)*jsz)+(ng-1)*jgz-distz <nz))	
		{ printf("geophones exceeds the computing zone!\n"); exit(1);}
	}
	else{
		if (!(gxbeg>=0 && gzbeg>=0 && gxbeg+(ng-1)*jgx<nx && gzbeg+(ng-1)*jgz<nz))	
		{ printf("geophones exceeds the computing zone!\n"); exit(1);}
	}
	cuda_set_sg<<<(ng+511)/512,512>>>(d_gxz, gxbeg, gzbeg, jgx, jgz, ng, nz);
	cudaMemset(d_bndr, 0, nt*2*(nz+nx)*sizeof(float));
	cudaMemset(d_dobs, 0, ng*nt*sizeof(float));
	cudaMemset(d_dcal, 0, ng*sizeof(float));
	cudaMemset(d_derr, 0, ns*ng*nt*sizeof(float));
	cudaMemset(d_g0, 0, nz*nx*sizeof(float));
	cudaMemset(d_g1, 0, nz*nx*sizeof(float));
	cudaMemset(d_cg, 0, nz*nx*sizeof(float));
	cudaMemset(d_lap, 0, nz*nx*sizeof(float));
	cudaMemset(d_sillum, 0, nz*nx*sizeof(float));
	cudaMemset(d_gillum, 0, nz*nx*sizeof(float));
	cudaMemset(d_pars, 0, 4*sizeof(float));
	cudaMemset(d_alpha1, 0, ng*sizeof(float));
	cudaMemset(d_alpha2, 0, ng*sizeof(float));
	cudaMemset(d_vtmp, 0, nz*nx*sizeof(float));
	cuda_init_bell<<<1,2*nbell+1>>>(d_bell);

	
	dim3 dimg=dim3((nz+Block_Size1-1)/Block_Size1, (nx+Block_Size2-1)/Block_Size2),dimb=dim3(Block_Size1, Block_Size2); 
	float dtx=dt/dx; 
	float dtz=dt/dz; 
	float _dz2=1.0/(dz*dz);
	float _dx2=1.0/(dx*dx);

	float *ptr=NULL;

	float mstimer = 0;// timer unit: millionseconds
	cudaEvent_t start, stop;
  	cudaEventCreate(&start);	
	cudaEventCreate(&stop);
	fp=fopen(shots_file,"rb");
	if (fp==NULL) { fprintf(stderr, "cannot open file %s:\n",__FILE__); exit(1);}
	for(int iter=0;iter<niter;iter++)
	{
		cudaEventRecord(start);
		rewind(fp);
		cudaMemcpy(d_g0, d_g1, nz*nx*sizeof(float), cudaMemcpyDeviceToDevice);
		cudaMemset(d_g1, 0, nz*nx*sizeof(float));
		cudaMemset(d_sillum, 0, nz*nx*sizeof(float));
		cudaMemset(d_gillum, 0, nz*nx*sizeof(float));
		for(int is=0;is<ns;is++)
		{
        		fread(dobs, sizeof(float), ng*nt, fp);
			matrix_transpose(dobs, nt, ng);
			cudaMemcpy(d_dobs, dobs, ng*nt*sizeof(float), cudaMemcpyHostToDevice);
			if (csdgather)	{
				gxbeg=sxbeg+is*jsx-distx;
				cuda_set_sg<<<(ng+511)/512, 512>>>(d_gxz, gxbeg, gzbeg, jgx, jgz, ng, nz);
			}
			wavefield_init(d_sp0, d_sp1, nz*nx);
			for(int it=0; it<nt; it++)
			{
				cuda_add_source<<<1,1>>>(d_sp1, &d_wlt[it], &d_sxz[is], 1, true);
				cuda_step_forward<<<dimg,dimb>>>(d_sp0, d_sp1, d_vv, dtz, dtx, nz, nx);
				ptr=d_sp0; d_sp0=d_sp1; d_sp1=ptr;

				cuda_record<<<(ng+511)/512, 512>>>(d_sp0, d_dcal, d_gxz, ng);
				cuda_cal_residuals<<<(ng+511)/512, 512>>>(d_dcal, &d_dobs[it*ng], &d_derr[is*ng*nt+it*ng], ng);
				cuda_rw_bndr<<<(2*nz+nx+512)/512,512>>>(&d_bndr[it*2*(nz+nx)], d_sp0, nz, nx, true);
			}
			//cudaMemcpy(dobs, d_dobs, ng*nt*sizeof(float), cudaMemcpyDeviceToHost);
			//matrix_transpose(dobs, ng, nt);
			//fwrite(dobs, sizeof(float), nt*ng, fp);

			ptr=d_sp0;d_sp0=d_sp1;d_sp1=ptr;
			wavefield_init(d_gp0, d_gp1, nz*nx);
			for(int it=nt-1; it>-1; it--)
			{
				cuda_rw_bndr<<<(2*nz+nx+255)/256,256>>>(&d_bndr[it*2*(nz+nx)], d_sp1, nz, nx, false);
				cuda_step_backward<<<dimg,dimb>>>(d_lap, d_sp0, d_sp1, d_vv, dtz, dtx, nz, nx);
				cuda_add_source<<<1,1>>>(d_sp1, &d_wlt[it], &d_sxz[is], 1, false);

				cuda_add_source<<<(ng+511)/512, 512>>>(d_gp1, &d_derr[is*ng*nt+it*ng], d_gxz, ng, true);
				cuda_step_forward<<<dimg,dimb>>>(d_gp0, d_gp1, d_vv, dtz, dtx, nz, nx);

				//cuda_cal_grad<<<dimg,dimb>>>(d_g1, d_sillum, d_gillum, d_lap, d_gp1, _dz2, _dx2, nz, nx);
				cuda_cal_gradient<<<dimg,dimb>>>(d_g1, d_sillum, d_lap, d_gp1, _dz2, _dx2, nz, nx);
				ptr=d_sp0; d_sp0=d_sp1; d_sp1=ptr;
				ptr=d_gp0; d_gp0=d_gp1; d_gp1=ptr;
			}
		}
		cuda_scale_gradient<<<dimg,dimb>>>(d_g1, d_vv, d_sillum, nz, nx);
		//cuda_scale_grad<<<dimg,dimb>>>(d_g1, d_vv, d_sillum, d_gillum, nz, nx);
		//cuda_gaussian_smoothz<<<dimg,dimb>>>(d_g1, d_sillum, d_bell, nz, nx);
		//cuda_gaussian_smoothx<<<dimg,dimb>>>(d_sillum, d_g1, d_bell, nz, nx);

		cuda_cal_objective<<<1, Block_Size>>>(&d_pars[0], d_derr, ns*ng*nt);
		report_par(&d_pars[0], "obj");

		FILE *fp1=fopen("grad.bin","wb");
		cudaMemcpy(v0, d_g1, nz*nx*sizeof(float), cudaMemcpyDeviceToHost);
		fwrite(v0, sizeof(float), nz*nx, fp1);
		fclose(fp1);

		if (iter>0) cuda_cal_beta<<<1, Block_Size>>>(&d_pars[1], d_g0, d_g1, d_cg, nz*nx); 
		cuda_cal_conjgrad<<<dimg, dimb>>>(d_g1, d_cg, &d_pars[1], nz, nx);
		cuda_cal_epsilon<<<1, Block_Size>>>(d_vv, d_cg, &d_pars[2], nz*nx);
		report_par(&d_pars[1], "beta");
		report_par(&d_pars[2], "epsil");

		rewind(fp);
		cudaMemset(d_alpha1, 0, ng*sizeof(float));
		cudaMemset(d_alpha2, 0, ng*sizeof(float));
		cuda_cal_vtmp<<<dimg, dimb>>>(d_vtmp, d_vv, d_cg, &d_pars[2], nz, nx);
		for(int is=0;is<ns;is++)
		{
        		fread(dobs, sizeof(float), ng*nt, fp);
			matrix_transpose(dobs, nt, ng);
			cudaMemcpy(d_dobs, dobs, ng*nt*sizeof(float), cudaMemcpyHostToDevice);
			if (csdgather)	{
				gxbeg=sxbeg+is*jsx-distx;
				cuda_set_sg<<<(ng+511)/512, 512>>>(d_gxz, gxbeg, gzbeg, jgx, jgz, ng, nz);
			}
			wavefield_init(d_sp0, d_sp1, nz*nx);
			for(int it=0; it<nt; it++)
			{
				cuda_add_source<<<1,1>>>(d_sp1, &d_wlt[it], &d_sxz[is], 1, true);
				cuda_step_forward<<<dimg,dimb>>>(d_sp0, d_sp1, d_vtmp, dtz, dtx, nz, nx);
				ptr=d_sp0; d_sp0=d_sp1; d_sp1=ptr;

				cuda_record<<<(ng+511)/512, 512>>>(d_sp0, d_dcal, d_gxz, ng);
				cuda_sum_alpha12<<<(ng+511)/512, 512>>>(d_alpha1, d_alpha2, d_dcal, &d_dobs[it*ng], &d_derr[is*ng*nt+it*ng], ng);
			}
		}
		cuda_cal_alpha<<<1,Block_Size>>>(&d_pars[3], d_alpha1, d_alpha2, &d_pars[2], ng);
		report_par(&d_pars[3], "alpha");
		cuda_update_vel<<<dimg,dimb>>>(d_vv, d_cg, &d_pars[3], nz, nx);


		FILE *fp2=fopen("newvel.bin","wb");
		cudaMemcpy(v0, d_vv, nz*nx*sizeof(float), cudaMemcpyDeviceToHost);
		fwrite(v0, sizeof(float), nz*nx, fp2);
		fclose(fp2);

		cudaEventRecord(stop);
  		cudaEventSynchronize(stop);
  		cudaEventElapsedTime(&mstimer, start, stop);
    		printf("iteration %d finished: %f (s)\n",iter+1, mstimer*1e-3);
	}
	fclose(fp);
	cudaEventDestroy(start);
	cudaEventDestroy(stop);

	free(v0);
	free(dobs);
	device_free();

	return 0;
}
