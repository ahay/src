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
#define EPS	1.0e-10
#endif

#define PI 	3.141592653589793
#define Block_Size1 16	/* 1st dim block size */
#define Block_Size2 16	/* 2nd dim block size */
#define Block_Size  512	/* vector computation blocklength */

#include "cuda_fwi_kernels.cu"

static bool csdgather;
static int niter,nz,nx, nz1,nx1,nt,ns,ng,sxbeg,szbeg,gxbeg,gzbeg,jsx,jsz,jgx,jgz;
static float dx, dz, fm, dt;

/* variables on host */
float 	*v0, *vv, *dobs;
/* variables on device */
int 	*d_sxz, *d_gxz;			
float 	*d_wlt, *d_vv, *d_illum, *d_lap, *d_vtmp, *d_sp0, *d_sp1, *d_gp0, *d_gp1,*d_bndr;
float	*d_dobs, *d_dcal, *d_derr, *d_g0, *d_g1, *d_cg, *d_pars, *d_alpha1, *d_alpha2;
/*
d_pars[0]: obj;
d_pars[1]: beta;
d_pars[2]: epsilon;
d_pars[3]: alpha;
d_alpha1[]: numerator of alpha, length=ng
d_alpha2[]: denominator of alpha, length=ng
*/

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

void matrix_transpose(float *matrix, int n1, int n2)
/*< matrix transpose >*/
{
	int i2, i1;
	float *tmp;
	if (!(tmp=(float*)malloc(n1*n2*sizeof(float))))
	 {sf_warning("out of memory!"); exit(1);}
	
	for(i2=0; i2<n2; i2++){
		for(i1=0; i1<n1; i1++){
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
	cudaMalloc(&d_bndr, nt*(2*nz+nx)*sizeof(float));
	cudaMalloc(&d_dobs, ng*nt*sizeof(float));
	cudaMalloc(&d_dcal, ng*sizeof(float));
	cudaMalloc(&d_derr, ns*ng*nt*sizeof(float));
	cudaMalloc(&d_g0, nz*nx*sizeof(float));
	cudaMalloc(&d_g1, nz*nx*sizeof(float));
	cudaMalloc(&d_cg, nz*nx*sizeof(float));
	cudaMalloc(&d_lap, nz*nx*sizeof(float));
	cudaMalloc(&d_illum, nz*nx*sizeof(float));
	cudaMalloc(&d_pars, 4*sizeof(float));
	cudaMalloc(&d_alpha1, ng*sizeof(float));
	cudaMalloc(&d_alpha2, ng*sizeof(float));
	cudaMalloc(&d_vtmp, nz*nx*sizeof(float));

    	cudaError_t err = cudaGetLastError ();
    	if (cudaSuccess != err) 
	sf_warning("Cuda error: Failed to allocate required memory!: %s", cudaGetErrorString(err));
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
	cudaFree(d_illum);
	cudaFree(d_pars);
	cudaFree(d_alpha1);
	cudaFree(d_alpha2);
	cudaFree(d_vtmp);

    	cudaError_t err = cudaGetLastError ();
    	if (cudaSuccess != err)
	sf_warning("Cuda error: Failed to free the allocated memory!: %s", cudaGetErrorString(err));
}

int main(int argc, char *argv[])
{
	bool verb;
	int is, it, iter, distx, distz, csd, rbell;
	float dtx, dtz, mstimer,amp, obj, beta, epsil, alpha;
	float *objval, *ptr=NULL;
	sf_file vinit, shots, vupdates, grads, objs, illums;

    	/* initialize Madagascar */
    	sf_init(argc,argv);

    	/*< set up I/O files >*/
    	vinit=sf_input ("in");   /* initial velocity model, unit=m/s */
	shots=sf_input("shots"); /* recorded shots from exact velocity model */
    	vupdates=sf_output("out"); /* updated velocity in iterations */ 
    	grads=sf_output("grads");  /* gradient in iterations */ 
	objs=sf_output("objs");/* values of objective function in iterations */
	illums=sf_output("illums");/* source illumination in iterations */

    	/* get parameters from velocity model and recorded shots */
	if (!sf_getbool("verb",&verb)) verb=true;
    	if (!sf_histint(vinit,"n1",&nz1)) sf_error("no n1");
    	if (!sf_histint(vinit,"n2",&nx1)) sf_error("no n2");
    	if (!sf_histfloat(vinit,"d1",&dz)) sf_error("no d1");
   	if (!sf_histfloat(vinit,"d2",&dx)) sf_error("no d2");

   	if (!sf_histint(shots,"n1",&nt)) sf_error("no nt");
	/* total modeling time steps */
   	if (!sf_histint(shots,"n2",&ng)) sf_error("no ng");
	/* total receivers in each shot */
   	if (!sf_histint(shots,"n3",&ns)) sf_error("no ns");
	/* number of shots */
   	if (!sf_histfloat(shots,"d1",&dt)) sf_error("no dt");
	/* time sampling interval */
   	if (!sf_histfloat(shots,"amp",&amp)) sf_error("no amp");
	/* maximum amplitude of ricker */
   	if (!sf_histfloat(shots,"fm",&fm)) sf_error("no fm");
	/* dominant freq of ricker */
   	if (!sf_histint(shots,"sxbeg",&sxbeg)) sf_error("no sxbeg");
	/* x-begining index of sources, starting from 0 */
   	if (!sf_histint(shots,"szbeg",&szbeg)) sf_error("no szbeg");
	/* x-begining index of sources, starting from 0 */
   	if (!sf_histint(shots,"gxbeg",&gxbeg)) sf_error("no gxbeg");
	/* x-begining index of receivers, starting from 0 */
   	if (!sf_histint(shots,"gzbeg",&gzbeg)) sf_error("no gzbeg");
	/* x-begining index of receivers, starting from 0 */
   	if (!sf_histint(shots,"jsx",&jsx)) sf_error("no jsx");
	/* source x-axis  jump interval  */
   	if (!sf_histint(shots,"jsz",&jsz)) sf_error("no jsz");
	/* source z-axis jump interval  */
   	if (!sf_histint(shots,"jgx",&jgx)) sf_error("no jgx");
	/* receiver x-axis jump interval  */
   	if (!sf_histint(shots,"jgz",&jgz)) sf_error("no jgz");
	/* receiver z-axis jump interval  */
   	if (!sf_histint(shots,"csdgather",&csd)) sf_error("csdgather or not required");
	/* default, common shot-gather; if n, record at every point*/
    	if (!sf_getint("niter",&niter))   niter=100;
	/* number of iterations */
	if (!sf_getint("rbell",&rbell))	  rbell=2;

	sf_putint(vupdates,"n1",nz1);	
	sf_putint(vupdates,"n2",nx1);
	sf_putfloat(vupdates,"d1",dz);
	sf_putfloat(vupdates,"d2",dx);
	sf_putstring(vupdates,"label1","Depth");
	sf_putstring(vupdates,"label2","Distance");
	sf_putstring(vupdates,"label3","Iteration");
	sf_putint(vupdates,"n3",niter);
	sf_putint(vupdates,"d3",1);
	sf_putint(vupdates,"o3",1);
	sf_putint(grads,"n1",nz1);	
	sf_putint(grads,"n2",nx1);
	sf_putint(grads,"n3",niter);
	sf_putfloat(grads,"d1",dz);
	sf_putfloat(grads,"d2",dx);
	sf_putint(grads,"d3",1);
	sf_putint(grads,"o3",1);
	sf_putstring(grads,"label1","Depth");
	sf_putstring(grads,"label2","Distance");
	sf_putstring(grads,"label3","Iteration");
	sf_putint(illums,"n1",nz1);	
	sf_putint(illums,"n2",nx1);
	sf_putfloat(illums,"d1",dz);
	sf_putfloat(illums,"d2",dx);
	sf_putint(illums,"n3",niter);
	sf_putint(illums,"d3",1);
	sf_putint(illums,"o3",1);
	sf_putint(objs,"n1",niter);
	sf_putint(objs,"n2",1);
	sf_putint(objs,"d1",1);
	sf_putint(objs,"o1",1);

	dtx=dt/dx; 
	dtz=dt/dz; 
	csdgather=(csd>0)?true:false;
	/* round the size up to multiples of Block size */
	nx=(int)((nx1+Block_Size1-1)/Block_Size1)*Block_Size1;
	nz=(int)((nz1+Block_Size2-1)/Block_Size2)*Block_Size2; 
	dim3 dimg=dim3(nz/Block_Size1, nx/Block_Size2), dimb=dim3(Block_Size1, Block_Size2); 

	v0=(float*)malloc(nz1*nx1*sizeof(float));
	vv=(float*)malloc(nz*nx*sizeof(float));
	dobs=(float*)malloc(ng*nt*sizeof(float));
	objval=(float*)malloc(niter*sizeof(float));
	sf_floatread(v0, nz1*nx1, vinit);
	expand(vv, v0, nz, nx, nz1, nx1);
	memset(dobs,0,ng*nt*sizeof(float));
	memset(objval,0,niter*sizeof(float));

    	cudaSetDevice(0);
    	cudaError_t err = cudaGetLastError();
    	if (cudaSuccess != err) 
	sf_warning("Cuda error: Failed to initialize device: %s\n", cudaGetErrorString(err));
	device_alloc(); 


	cudaMemcpy(d_vv, vv, nz*nx*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemset(d_sp0,0,nz*nx*sizeof(float));
	cudaMemset(d_sp1,0,nz*nx*sizeof(float));
	cudaMemset(d_gp0,0,nz*nx*sizeof(float));
	cudaMemset(d_gp1,0,nz*nx*sizeof(float));
	cuda_ricker_wavelet<<<(nt+511)/512,512>>>(d_wlt, amp, fm, dt, nt);
	if (!(sxbeg>=0 && szbeg>=0 && sxbeg+(ns-1)*jsx<nx1 && szbeg+(ns-1)*jsz<nz1))	
	{ sf_warning("sources exceeds the computing zone!\n"); exit(1);}

	cuda_set_sg<<<(ns+511)/512,512>>>(d_sxz, sxbeg, szbeg, jsx, jsz, ns, nz);
	distx=sxbeg-gxbeg;
	distz=szbeg-gzbeg;
	if (csdgather)	{
		if (!(gxbeg>=0 && gzbeg>=0 && gxbeg+(ng-1)*jgx<nx1 && gzbeg+(ng-1)*jgz<nz1 &&
		(sxbeg+(ns-1)*jsx)+(ng-1)*jgx-distx <nx1  && (szbeg+(ns-1)*jsz)+(ng-1)*jgz-distz <nz1))	
		{ sf_warning("geophones exceeds the computing zone!\n"); exit(1);}
	}
	else{
		if (!(gxbeg>=0 && gzbeg>=0 && gxbeg+(ng-1)*jgx<nx1 && gzbeg+(ng-1)*jgz<nz1))	
		{ sf_warning("geophones exceeds the computing zone!\n"); exit(1);}
	}
	cuda_set_sg<<<(ng+511)/512,512>>>(d_gxz, gxbeg, gzbeg, jgx, jgz, ng, nz);
	cudaMemset(d_bndr, 0, nt*(2*nz+nx)*sizeof(float));
	cudaMemset(d_dobs, 0, ng*nt*sizeof(float));
	cudaMemset(d_dcal, 0, ng*sizeof(float));
	cudaMemset(d_derr, 0, ns*ng*nt*sizeof(float));
	cudaMemset(d_g0, 0, nz*nx*sizeof(float));
	cudaMemset(d_g1, 0, nz*nx*sizeof(float));
	cudaMemset(d_cg, 0, nz*nx*sizeof(float));
	cudaMemset(d_lap, 0, nz*nx*sizeof(float));
	cudaMemset(d_illum, 0, nz*nx*sizeof(float));
	cudaMemset(d_pars, 0, 4*sizeof(float));
	cudaMemset(d_alpha1, 0, ng*sizeof(float));
	cudaMemset(d_alpha2, 0, ng*sizeof(float));
	cudaMemset(d_vtmp, 0, nz*nx*sizeof(float));
	
	cudaEvent_t start, stop;
  	cudaEventCreate(&start);	
	cudaEventCreate(&stop);
	for(iter=0; iter<niter; iter++)
	{
		cudaEventRecord(start);

		sf_seek(shots, 0L, SEEK_SET);
		cudaMemcpy(d_g0, d_g1, nz*nx*sizeof(float), cudaMemcpyDeviceToDevice);
		cudaMemset(d_g1, 0, nz*nx*sizeof(float));
		cudaMemset(d_illum, 0, nz*nx*sizeof(float));
		for(is=0;is<ns;is++)
		{
			sf_floatread(dobs, ng*nt, shots);
			matrix_transpose(dobs, nt, ng);
			cudaMemcpy(d_dobs, dobs, ng*nt*sizeof(float), cudaMemcpyHostToDevice);
			if (csdgather)	{
				gxbeg=sxbeg+is*jsx-distx;
				cuda_set_sg<<<(ng+511)/512, 512>>>(d_gxz, gxbeg, gzbeg, jgx, jgz, ng, nz);
			}
			cudaMemset(d_sp0,0,nz*nx*sizeof(float));
			cudaMemset(d_sp1,0,nz*nx*sizeof(float));
			for(it=0; it<nt; it++)
			{
				cuda_add_source<<<1,1>>>(d_sp1, &d_wlt[it], &d_sxz[is], 1, true);
				cuda_step_forward<<<dimg,dimb>>>(d_sp0, d_sp1, d_vv, dtz, dtx, nz, nx);
				ptr=d_sp0; d_sp0=d_sp1; d_sp1=ptr;

				cuda_record<<<(ng+511)/512, 512>>>(d_sp0, d_dcal, d_gxz, ng);
				cuda_cal_residuals<<<(ng+511)/512, 512>>>(d_dcal, &d_dobs[it*ng], &d_derr[is*ng*nt+it*ng], ng);
				cuda_rw_bndr<<<(2*nz+nx+511)/512,512>>>(&d_bndr[it*(2*nz+nx)], d_sp0, nz, nx, true);
			}

			ptr=d_sp0;d_sp0=d_sp1;d_sp1=ptr;
			cudaMemset(d_gp0,0,nz*nx*sizeof(float));
			cudaMemset(d_gp1,0,nz*nx*sizeof(float));
			for(it=nt-1; it>-1; it--)
			{
				cuda_rw_bndr<<<(2*nz+nx+255)/256,256>>>(&d_bndr[it*(2*nz+nx)], d_sp1, nz, nx, false);
				cuda_step_backward<<<dimg,dimb>>>(d_lap, d_sp0, d_sp1, d_vv, dtz, dtx, nz, nx);
				cuda_add_source<<<1,1>>>(d_sp1, &d_wlt[it], &d_sxz[is], 1, false);

				cuda_add_source<<<(ng+511)/512, 512>>>(d_gp1, &d_derr[is*ng*nt+it*ng], d_gxz, ng, true);
				cuda_step_forward<<<dimg,dimb>>>(d_gp0, d_gp1, d_vv, dtz, dtx, nz, nx);

				cuda_cal_gradient<<<dimg,dimb>>>(d_g1, d_illum, d_lap, d_gp1, nz, nx);
				ptr=d_sp0; d_sp0=d_sp1; d_sp1=ptr;
				ptr=d_gp0; d_gp0=d_gp1; d_gp1=ptr;
			}
		}
		cuda_cal_objective<<<1, Block_Size>>>(&d_pars[0], d_derr, ns*ng*nt);
		cudaMemcpy(&obj, &d_pars[0], sizeof(float), cudaMemcpyDeviceToHost);

		cudaMemcpy(vv, d_illum, nz*nx*sizeof(float), cudaMemcpyDeviceToHost);
		window(v0, vv, nz, nx, nz1, nx1);
		sf_floatwrite(v0, nz1*nx1, illums);

		cuda_scale_gradient<<<dimg,dimb>>>(d_g1, d_vv, d_illum, nz, nx);
		cudaMemcpy(vv, d_g1, nz*nx*sizeof(float), cudaMemcpyDeviceToHost);
		window(v0, vv, nz, nx, nz1, nx1);
		sf_floatwrite(v0, nz1*nx1, grads);
		cuda_bell_smoothz<<<dimg,dimb>>>(d_g1, d_illum, rbell, nz, nx);
		cuda_bell_smoothx<<<dimg,dimb>>>(d_illum, d_g1, rbell, nz, nx);

		if (iter>0) cuda_cal_beta<<<1, Block_Size>>>(&d_pars[1], d_g0, d_g1, d_cg, nz*nx); 
		cudaMemcpy(&beta, &d_pars[1], sizeof(float), cudaMemcpyDeviceToHost);
		cuda_cal_conjgrad<<<dimg, dimb>>>(d_g1, d_cg, beta, nz, nx);

		cuda_cal_epsilon<<<1, Block_Size>>>(d_vv, d_cg, &d_pars[2], nz*nx);
		cudaMemcpy(&epsil, &d_pars[2], sizeof(float), cudaMemcpyDeviceToHost);

		sf_seek(shots, 0L, SEEK_SET);
		cudaMemset(d_alpha1, 0, ng*sizeof(float));
		cudaMemset(d_alpha2, 0, ng*sizeof(float));
		cuda_cal_vtmp<<<dimg, dimb>>>(d_vtmp, d_vv, d_cg, epsil, nz, nx);
		for(is=0;is<ns;is++)
		{
			sf_floatread(dobs, ng*nt, shots);
			matrix_transpose(dobs, nt, ng);
			cudaMemcpy(d_dobs, dobs, ng*nt*sizeof(float), cudaMemcpyHostToDevice);
			if (csdgather)	{
				gxbeg=sxbeg+is*jsx-distx;
				cuda_set_sg<<<(ng+511)/512, 512>>>(d_gxz, gxbeg, gzbeg, jgx, jgz, ng, nz);
			}
			cudaMemset(d_sp0,0,nz*nx*sizeof(float));
			cudaMemset(d_sp1,0,nz*nx*sizeof(float));
			for(it=0; it<nt; it++)
			{
				cuda_add_source<<<1,1>>>(d_sp1, &d_wlt[it], &d_sxz[is], 1, true);
				cuda_step_forward<<<dimg,dimb>>>(d_sp0, d_sp1, d_vtmp, dtz, dtx, nz, nx);
				ptr=d_sp0; d_sp0=d_sp1; d_sp1=ptr;

				cuda_record<<<(ng+511)/512, 512>>>(d_sp0, d_dcal, d_gxz, ng);
				cuda_sum_alpha12<<<(ng+511)/512, 512>>>(d_alpha1, d_alpha2, d_dcal, &d_dobs[it*ng], &d_derr[is*ng*nt+it*ng], ng);
			}
		}
		cuda_cal_alpha<<<1,Block_Size>>>(&d_pars[3], d_alpha1, d_alpha2, epsil, ng);
		cudaMemcpy(&alpha, &d_pars[3], sizeof(float), cudaMemcpyDeviceToHost);

		cuda_update_vel<<<dimg,dimb>>>(d_vv, d_cg, alpha, nz, nx);
		cudaMemcpy(vv, d_vv, nz*nx*sizeof(float), cudaMemcpyDeviceToHost);
		window(v0, vv, nz, nx, nz1, nx1);
		sf_floatwrite(v0, nz1*nx1, vupdates);

		cudaEventRecord(stop);
  		cudaEventSynchronize(stop);
  		cudaEventElapsedTime(&mstimer, start, stop);

		if(verb) {
			objval[iter]=obj;
			sf_warning("obj=%f  beta=%f  epsil=%f  alpha=%f",obj, beta, epsil, alpha);
			sf_warning("iteration %d finished: %f (s)",iter+1, mstimer*1e-3);
		}
	}
	cudaEventDestroy(start);
	cudaEventDestroy(stop);

	sf_floatwrite(objval,iter,objs);
	sf_fileclose(shots);

	free(v0);
	free(vv);
	free(dobs);
	free(objval);
	device_free();

	return 0;
}
