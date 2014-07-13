/* Forward-backword exact reconstruction using boundary saving

Note: 	It is used as a demonstration that we can reconstruct the modeled
	wavefield exactly via boundary saving.
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

#define PI 	SF_PI
#define Block_Size1 16		
#define Block_Size2 16		
#define Block_Size  512		
#define nbell	2		

#include "cuda_fwi_kernels.cu"


void sf_check_gpu_error (const char *msg) 
/*< check GPU errors >*/
{
    cudaError_t err = cudaGetLastError ();
    if (cudaSuccess != err) { 
	sf_error ("Cuda error: %s: %s", msg, cudaGetErrorString (err)); 
	exit(0);   
    }
}

void matrix_transpose(float *matrix, int n1, int n2)
/*< matrix transpose >*/
{
	int i1, i2;
	float *tmp=(float*)malloc(n1*n2*sizeof(float));
	if (tmp==NULL) {printf("out of memory!"); exit(1);}
	for(i2=0; i2<n2; i2++){
		for(i1=0; i1<n1; i1++){
			tmp[i2+n2*i1]=matrix[i1+n1*i2];
		}
	}
	memcpy(matrix, tmp, n1*n2*sizeof(float));
	free(tmp);
}

void expand(float*vv, float *v0, int nz, int nx, int nz1, int nx1)
/*<expand model from nz1*nx1 to nz*nx >*/
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
/*<window a model from nz*nx to nz1*nx1 >*/
{
	int i1, i2;

	for(i2=0; i2<nx1; i2++)
	for(i1=0; i1<nz1; i1++)
		  v0[i1+i2*nz1]=vv[i1+nz*i2];
}

int main(int argc, char *argv[])
{
	/* variables on host */
	bool csdgather;
	int is, ft, jt, it, distx, distz;
	int nz,nx,nz1,nx1,nt,ns,ng,sxbeg,szbeg,gxbeg,gzbeg,jsx,jsz,jgx,jgz;
	float dx, dz, fm, dt, dtx,dtz,amp;
	float 	*v0, *dobs, *vv, *ptr=NULL;
	sf_file vinit, Fw1, Fw2;
	/* variables on device */
	int 	*d_sxz, *d_gxz;			
	float 	*d_wlt, *d_vv, *d_sp0, *d_sp1, *d_lap, *d_illum, *d_dobs, *d_bndr;

    	/* initialize Madagascar */
    	sf_init(argc,argv);

    	/*< set up I/O files >*/
    	vinit=sf_input ("in");/* initial velocity model, unit=m/s */
	Fw1 = sf_output("out");/* forward wavefield snaps */
	Fw2 = sf_output("back");/* backward wavefield snaps */

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
    	if (!sf_getint("ns",&ns))   ns=1;	
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
	if (!sf_getbool("csdgather",&csdgather)) csdgather=true;
	/* default, common shot-gather; if n, record at every point*/
   	if (!sf_getint("ft",&ft)) ft=0; 
	/* first recorded time */
    	if (!sf_getint("jt",&jt)) jt=1;	
	/* time interval */

	sf_putint(Fw1,"n1",nz1);
	sf_putint(Fw1,"n2",nx1);
    	sf_putint(Fw1,"n3",(nt-ft)/jt);
    	sf_putfloat(Fw1,"d3",jt*dt);
    	sf_putfloat(Fw1,"o3",ft*dt);
	sf_putint(Fw2,"n1",nz1);
	sf_putint(Fw2,"n2",nx1);
    	sf_putint(Fw2,"n3",(nt-ft)/jt);
    	sf_putfloat(Fw2,"d3",-jt*dt);
    	sf_putfloat(Fw2,"o3",nt*dt);

	dtx=dt/dx; 
	dtz=dt/dz; 
	/* round the size up to multiples of Block size */
	nx=(int)((nx1+Block_Size1-1)/Block_Size1)*Block_Size1;
	nz=(int)((nz1+Block_Size2-1)/Block_Size2)*Block_Size2;

	v0=(float*)malloc(nz1*nx1*sizeof(float));
	vv=(float*)malloc(nz*nx*sizeof(float));
	dobs=(float*)malloc(ng*nt*sizeof(float));
	sf_floatread(v0,nz1*nx1,vinit);
	expand(vv, v0, nz, nx, nz1, nx1);
	memset(dobs,0,ng*nt*sizeof(float));

    	cudaSetDevice(0);
	sf_check_gpu_error("Failed to initialize device!");
	cudaMalloc(&d_vv, nz*nx*sizeof(float));
	cudaMalloc(&d_sp0, nz*nx*sizeof(float));
	cudaMalloc(&d_sp1, nz*nx*sizeof(float));
	cudaMalloc(&d_lap, nz*nx*sizeof(float));
	cudaMalloc(&d_illum, nz*nx*sizeof(float));
	cudaMalloc(&d_wlt, nt*sizeof(float));
	cudaMalloc(&d_sxz, ns*sizeof(float));
	cudaMalloc(&d_gxz, ng*sizeof(float));
	cudaMalloc(&d_dobs, ng*nt*sizeof(float));
	cudaMalloc(&d_bndr, nt*(2*nx+nz)*sizeof(float));
	sf_check_gpu_error("Failed to allocate required memory on device!");

	dim3 dimg=dim3(nz/Block_Size1, nx/Block_Size2),dimb=dim3(Block_Size1, Block_Size2); 

	cudaMemcpy(d_vv, vv, nz*nx*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemset(d_sp0,0,nz*nx*sizeof(float));
	cudaMemset(d_sp1,0,nz*nx*sizeof(float));
	cuda_ricker_wavelet<<<(nt+511)/512,512>>>(d_wlt,amp, fm, dt, nt);
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
	
	for(is=0; is<ns; is++)
	{
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
			cuda_rw_bndr<<<((2*nz+nx)+511)/512,512>>>(&d_bndr[it*(2*nz+nx)], d_sp0, nz, nx, true);

			if(it>=ft)
			{
				cudaMemcpy(vv, d_sp0, nz*nx*sizeof(float), cudaMemcpyDeviceToHost);
				window(v0, vv, nz, nx, nz1, nx1);
				sf_floatwrite(v0,nz1*nx1,Fw1);
			}
		}

		ptr=d_sp0;d_sp0=d_sp1;d_sp1=ptr;
		for(it=nt-1; it>-1; it--)
		{
			if(it>=ft)
			{
				cudaMemcpy(vv, d_sp1, nz*nx*sizeof(float), cudaMemcpyDeviceToHost);
				window(v0, vv, nz, nx, nz1, nx1);
				sf_floatwrite(v0,nz1*nx1,Fw2);
			}

			cuda_rw_bndr<<<(2*(nz+nx)+511)/512,512>>>(&d_bndr[it*(2*nz+nx)], d_sp1, nz, nx, false);
			cuda_step_backward<<<dimg,dimb>>>(d_illum, d_lap, d_sp0, d_sp1, d_vv, dtz, dtx, nz, nx);
			cuda_add_source<<<1,1>>>(d_sp1, &d_wlt[it], &d_sxz[is], 1, false);

			ptr=d_sp0; d_sp0=d_sp1; d_sp1=ptr;
		}
	}

	/* free memory on device */
	cudaFree(d_vv);
	cudaFree(d_sp0);
	cudaFree(d_sp1);
	cudaFree(d_lap);
	cudaFree(d_illum);
	cudaFree(d_wlt);
	cudaFree(d_sxz);
	cudaFree(d_gxz);
	cudaFree(d_dobs);
	cudaFree(d_bndr);
	sf_check_gpu_error("Failed to free the allocated memory!");
	/* free memory on host */
	free(v0);
	free(vv);
	free(dobs);

	return 0;
}

