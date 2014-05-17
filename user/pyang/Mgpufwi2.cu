/* CUDA based FWI using PML absorbing boundary condition

Note: 	You can try other complex boundary condition but we do not
	recommend to do so. The main reason is that FWI is to recover
	the low-frequency information of the earth model. Low-freq 
	means that exact absorbing is not necessarilly needed. The 
	result will be improved with the optimization precedure. 
	Furthermore, complex boundary condition (such as sponge ABC or
	PML) implies more computational cost, which will degrade the
	efficiency of FWI. 
   	coordinate configuration of dobsmic data:
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

	 2nd-order stability condition:	min(dx, dz)>sqrt(2)*dt*max(v)
	 numerical dispersion condition:	max(dx, dz)<min(v)/(10*fmax)
*/
/*
  Copyright (C) 2013  Xi'an Jiaotong University (Pengliang Yang)
    Email: ypl.2100@gmail.com	
    Acknowledgement: This code is written with the help of Baoli Wang.

References:
    [1] Tarantola, Albert. "Inversion of seismic reflection data in the 
	acoustic approximation." Geophysics 49.8 (1984): 1259-1266.
    [2] Pica, A., J. P. Diet, and A. Tarantola. "Nonlinear inversion 
	of seismic reflection data in a laterally invariant medium." 
	Geophysics 55.3 (1990): 284-292.
    [3] Hager, William W., and Hongchao Zhang. "A survey of nonlinear
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
#define EPS	1.0e-15f
#endif

#define PI 	3.141592653589793f
#define Block_Size1 16		// 1st dim block size
#define Block_Size2 16		// 2nd dim block size
#define Block_Size  512
static const int npml=32;

#include "cuda_fwi2_kernels.cu"

static bool 	csdgather;
static int 	NJ,niter,nz,nx,nnx,nnz,nz1,nx1,nt,nt_h,ns,ng,sxbeg,szbeg,gxbeg,gzbeg,jsx,jsz,jgx,jgz;
static float 	dx, dz, _dx,_dz, fm, dt;
static dim3 	dimg, dimb, dimgx1, dimbx1, dimgx2, dimbx2, dimgz1, dimbz1, dimgz2, dimbz2;

// variables on host
float 	*dobs, *v0, *vel, *p;
// variables on device
int 	*d_sxz, *d_gxz;
float 	*d_sp0, *d_sp1, *d_svx, *d_svz, *d_gp0, *d_gp1, *d_gvx, *d_gvz;
float 	*d_bx1, *d_bx2, *d_bz1, *d_bz2, *d_convpx, *d_convpz, *d_convvx, *d_convvz;
float	*d_wlt, *d_vv, *d_vtmp, *h_boundary, *d_boundary, *d_pars;			
float 	*d_dcal, *d_dobs, *d_derr, *d_g0, *d_g1, *d_cg, *d_power, *d_alpha1, *d_alpha2;	
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
	if (tmp==NULL) {sf_warning("out of memory!"); exit(1);}
	for(int i2=0; i2<n2; i2++){
		for(int i1=0; i1<n1; i1++){
			tmp[i2+n2*i1]=matrix[i1+n1*i2];
		}
	}
	memcpy(matrix, tmp, n1*n2*sizeof(float));
	free(tmp);
}

void device_alloc()
{
    	cudaMalloc(&d_sxz,	ns*sizeof(int));
    	cudaMalloc(&d_gxz, 	ng*sizeof(int));
	cudaMalloc(&d_wlt,	nt*sizeof(float));
    	cudaMalloc(&d_vv, 	nnz*nnx*sizeof(float));
    	cudaMalloc(&d_sp0, 	nnz*nnx*sizeof(float));
    	cudaMalloc(&d_sp1, 	nnz*nnx*sizeof(float));
    	cudaMalloc(&d_svx, 	nnz*nnx*sizeof(float));
    	cudaMalloc(&d_svz, 	nnz*nnx*sizeof(float));
    	cudaMalloc(&d_gp0, 	nnz*nnx*sizeof(float));
    	cudaMalloc(&d_gp1, 	nnz*nnx*sizeof(float));
    	cudaMalloc(&d_gvx, 	nnz*nnx*sizeof(float));
    	cudaMalloc(&d_gvz, 	nnz*nnx*sizeof(float));
    	cudaMalloc(&d_bx1, 	2*npml*nnz*sizeof(float));
    	cudaMalloc(&d_bz1, 	2*npml*nnx*sizeof(float));
    	cudaMalloc(&d_bx2, 	2*npml*nnz*sizeof(float));
    	cudaMalloc(&d_bz2, 	2*npml*nnx*sizeof(float));
    	cudaMalloc(&d_convvx, 	2*npml*nnz*sizeof(float));
    	cudaMalloc(&d_convvz, 	2*npml*nnx*sizeof(float));
    	cudaMalloc(&d_convpx, 	2*npml*nnz*sizeof(float));
    	cudaMalloc(&d_convpz, 	2*npml*nnx*sizeof(float));
	cudaHostAlloc(&h_boundary, nt_h*2*(NJ-1)*(nx+nz)*sizeof(float), cudaHostAllocMapped);	
	cudaMalloc(&d_boundary, (nt-nt_h)*2*(NJ-1)*(nx+nz)*sizeof(float));
    	cudaMalloc(&d_dcal, 	ng*sizeof(float));
    	cudaMalloc(&d_dobs, 	ng*nt*sizeof(float));
    	cudaMalloc(&d_derr, 	ns*ng*nt*sizeof(float));
    	cudaMalloc(&d_g0, 	nnz*nnx*sizeof(float));
    	cudaMalloc(&d_g1, 	nnz*nnx*sizeof(float));
    	cudaMalloc(&d_cg, 	nnz*nnx*sizeof(float));
	cudaMalloc(&d_power, 	nnz*nnx*sizeof(float));
	cudaMalloc(&d_vtmp, 	nnz*nnx*sizeof(float));
	cudaMalloc(&d_alpha1, 	ng*sizeof(float));
	cudaMalloc(&d_alpha2, 	ng*sizeof(float));
	cudaMalloc(&d_pars, 	4*sizeof(float));
}

void device_free()
{
    	cudaFree(d_sxz);
    	cudaFree(d_gxz);
	cudaFree(d_wlt);
    	cudaFree(d_vv);
    	cudaFree(d_sp0);
    	cudaFree(d_sp1);
    	cudaFree(d_svx);
    	cudaFree(d_svz);
    	cudaFree(d_gp0);
    	cudaFree(d_gp1);
    	cudaFree(d_gvx);
    	cudaFree(d_gvz);
	cudaFree(d_bx1);
	cudaFree(d_bx2);
	cudaFree(d_bz1);
	cudaFree(d_bz2);
    	cudaFree(d_convvx);
    	cudaFree(d_convvz);
    	cudaFree(d_convpx);
    	cudaFree(d_convpz);
	cudaFreeHost(h_boundary);
    	cudaFree(d_boundary);
    	cudaFree(d_dcal);
    	cudaFree(d_dobs);
    	cudaFree(d_derr);
    	cudaFree(d_g0);
	cudaFree(d_g1);
	cudaFree(d_cg);
	cudaFree(d_power);
	cudaFree(d_vtmp);
    	cudaFree(d_alpha1);
    	cudaFree(d_alpha2);
    	cudaFree(d_pars);
}

void wavefield_init(float *d_p0, float *d_p1, float *d_vx, float *d_vz, float *d_convpx, float *d_convpz, float *d_convvx, float *d_convvz)
{
	cudaMemset(d_p0, 	0,	nnz*nnx*sizeof(float));
	cudaMemset(d_p1, 	0,	nnz*nnx*sizeof(float));
	cudaMemset(d_vx, 	0,	nnz*nnx*sizeof(float));
	cudaMemset(d_vz, 	0,	nnz*nnx*sizeof(float));
	cudaMemset(d_convpx, 	0,	2*npml*nnz*sizeof(float));
	cudaMemset(d_convpz, 	0,	2*npml*nnx*sizeof(float));
	cudaMemset(d_convvx, 	0,	2*npml*nnz*sizeof(float));
	cudaMemset(d_convvz, 	0,	2*npml*nnx*sizeof(float));
}

void step_forward(float *d_vv, float *d_p0, float *d_p1, float *d_vx, float *d_vz, float *d_convvx, 
	float *d_convvz, float *d_convpx, float *d_convpz, float *d_bx1, float *d_bz1, float *d_bx2, float *d_bz2)
{
	// p0: p{it-1}; 
	// p1: p{it}; 
	// p{it+1}-->p0
	if (NJ==2)
	{
		cuda_forward_v_2<<<dimg, dimb>>>(d_p1, d_vx, d_vz, _dx, _dz, npml, nnz, nnx);
		cuda_PML_vz_2<<<dimgz1, dimbz1>>>(d_p1, d_convpz, d_bz2, d_vz, _dz, npml, nnz, nnx);
		cuda_PML_vx_2<<<dimgx1, dimbx1>>>(d_p1, d_convpx, d_bx2, d_vx, _dx, npml, nnz, nnx);
		cuda_forward_p_2<<<dimg, dimb>>>(d_vv, d_p0, d_p1, d_vx, d_vz, dt, _dx, _dz, npml, nnz, nnx);
		cuda_PML_pz_2<<<dimgz1, dimbz1>>>(d_vv, d_p0, d_convvz, d_bz1, d_vz, dt, _dz, npml, nnz, nnx);
		cuda_PML_px_2<<<dimgx1, dimbx1>>>(d_vv, d_p0, d_convvx, d_bx1, d_vx, dt, _dx, npml, nnz, nnx);
	}
	else if (NJ==4)
	{
		cuda_forward_v_4<<<dimg, dimb>>>(d_p1, d_vx, d_vz, _dx, _dz, npml, nnz, nnx);
		cuda_PML_vz_4<<<dimgz1, dimbz1>>>(d_p1, d_convpz, d_bz2, d_vz, _dz, npml, nnz, nnx);
		cuda_PML_vx_4<<<dimgx1, dimbx1>>>(d_p1, d_convpx, d_bx2, d_vx, _dx, npml, nnz, nnx);
		cuda_forward_p_4<<<dimg, dimb>>>(d_vv, d_p0, d_p1, d_vx, d_vz, dt, _dx, _dz, npml, nnz, nnx);
		cuda_PML_pz_4<<<dimgz1, dimbz1>>>(d_vv, d_p0, d_convvz, d_bz1, d_vz, dt, _dz, npml, nnz, nnx);
		cuda_PML_px_4<<<dimgx1, dimbx1>>>(d_vv, d_p0, d_convvx, d_bx1, d_vx, dt, _dx, npml, nnz, nnx);
	}
	else if (NJ==6)
	{
		cuda_forward_v_6<<<dimg, dimb>>>(d_p1, d_vx, d_vz, _dx, _dz, npml, nnz, nnx);
		cuda_PML_vz_6<<<dimgz1, dimbz1>>>(d_p1, d_convpz, d_bz2, d_vz, _dz, npml, nnz, nnx);
		cuda_PML_vx_6<<<dimgx1, dimbx1>>>(d_p1, d_convpx, d_bx2, d_vx, _dx, npml, nnz, nnx);
		cuda_forward_p_6<<<dimg, dimb>>>(d_vv, d_p0, d_p1, d_vx, d_vz, dt, _dx, _dz, npml, nnz, nnx);
		cuda_PML_pz_6<<<dimgz1, dimbz1>>>(d_vv, d_p0, d_convvz, d_bz1, d_vz, dt, _dz, npml, nnz, nnx);
		cuda_PML_px_6<<<dimgx1, dimbx1>>>(d_vv, d_p0, d_convvx, d_bx1, d_vx, dt, _dx, npml, nnz, nnx);
	}
	else if (NJ==8)
	{
		cuda_forward_v_8<<<dimg, dimb>>>(d_p1, d_vx, d_vz, _dx, _dz, npml, nnz, nnx);
		cuda_PML_vz_8<<<dimgz1, dimbz1>>>(d_p1, d_convpz, d_bz2, d_vz, _dz, npml, nnz, nnx);
		cuda_PML_vx_8<<<dimgx1, dimbx1>>>(d_p1, d_convpx, d_bx2, d_vx, _dx, npml, nnz, nnx);
		cuda_forward_p_8<<<dimg, dimb>>>(d_vv, d_p0, d_p1, d_vx, d_vz, dt, _dx, _dz, npml, nnz, nnx);
		cuda_PML_pz_8<<<dimgz1, dimbz1>>>(d_vv, d_p0, d_convvz, d_bz1, d_vz, dt, _dz, npml, nnz, nnx);
		cuda_PML_px_8<<<dimgx1, dimbx1>>>(d_vv, d_p0, d_convvx, d_bx1, d_vx, dt, _dx, npml, nnz, nnx);
	}
	else if (NJ==10)
	{
		cuda_forward_v_10<<<dimg, dimb>>>(d_p1, d_vx, d_vz, _dx, _dz, npml, nnz, nnx);
		cuda_PML_vz_10<<<dimgz1, dimbz1>>>(d_p1, d_convpz, d_bz2, d_vz, _dz, npml, nnz, nnx);
		cuda_PML_vx_10<<<dimgx1, dimbx1>>>(d_p1, d_convpx, d_bx2, d_vx, _dx, npml, nnz, nnx);
		cuda_forward_p_10<<<dimg, dimb>>>(d_vv, d_p0, d_p1, d_vx, d_vz, dt, _dx, _dz, npml, nnz, nnx);
		cuda_PML_pz_10<<<dimgz1, dimbz1>>>(d_vv, d_p0, d_convvz, d_bz1, d_vz, dt, _dz, npml, nnz, nnx);
		cuda_PML_px_10<<<dimgx1, dimbx1>>>(d_vv, d_p0, d_convvx, d_bx1, d_vx, dt, _dx, npml, nnz, nnx);
	}
}

void step_backward(float *d_vv, float *d_p0, float *d_p1, float *d_vx, float *d_vz)
{
	// p0: p{it-1}; 
	// p1: p{it}; 
	// p{it+1}-->p0
	if (NJ==2)
	{
		cuda_forward_v_2<<<dimg, dimb>>>(d_p1, d_vx, d_vz, _dx, _dz, npml, nnz, nnx);
		cuda_forward_p_2<<<dimg, dimb>>>(d_vv, d_p0, d_p1, d_vx, d_vz, dt, _dx, _dz, npml, nnz, nnx);
	}
	else if (NJ==4)
	{
		cuda_forward_v_4<<<dimg, dimb>>>(d_p1, d_vx, d_vz, _dx, _dz, npml, nnz, nnx);
		cuda_forward_p_4<<<dimg, dimb>>>(d_vv, d_p0, d_p1, d_vx, d_vz, dt, _dx, _dz, npml, nnz, nnx);
	}
	else if (NJ==6)
	{
		cuda_forward_v_6<<<dimg, dimb>>>(d_p1, d_vx, d_vz, _dx, _dz, npml, nnz, nnx);
		cuda_forward_p_6<<<dimg, dimb>>>(d_vv, d_p0, d_p1, d_vx, d_vz, dt, _dx, _dz, npml, nnz, nnx);
	}
	else if (NJ==8)
	{
		cuda_forward_v_8<<<dimg, dimb>>>(d_p1, d_vx, d_vz, _dx, _dz, npml, nnz, nnx);
		cuda_forward_p_8<<<dimg, dimb>>>(d_vv, d_p0, d_p1, d_vx, d_vz, dt, _dx, _dz, npml, nnz, nnx);
	}
	else if (NJ==10)
	{
		cuda_forward_v_10<<<dimg, dimb>>>(d_p1, d_vx, d_vz, _dx, _dz, npml, nnz, nnx);
		cuda_forward_p_10<<<dimg, dimb>>>(d_vv, d_p0, d_p1, d_vx, d_vz, dt, _dx, _dz, npml, nnz, nnx);
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
	if ( 	((NJ==2) &&(fm>=minvel/(10*MAX(dx,dz))))||
		((NJ==4) &&(fm>=minvel/(5*MAX(dx,dz))))	)
	sf_warning("Non-dispersion relation not satisfied!");
}

// a: size=nz1*nx1; b: size=nnz*nnx; 
void expand(float *b, float *a, int npml, int nnz, int nnx, int nz1, int nx1)
{/*< expand domain of 'a' to 'b' >*/
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

// from: size=nnz*nnx; to: size=nz1*nx1;
void window(float *to, float *from, int npml, int nnz, int nnx, int nz1, int nx1)
{
	int ix, iz;
	for(ix=0;ix<nx1;ix++){
		for(iz=0;iz<nz1;iz++){
			to[iz+ix*nz1]=from[(iz+npml)+(ix+npml)*nnz];
		}
	}
}



int main(int argc, char *argv[])
{
	bool verb;
	int is, kt, iter, distx, distz, csd;
	float phost, mstimer,amp, obj, beta, epsil, alpha;
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
    	if (!sf_getint("order",&NJ))   NJ=4;
	/* order of finite difference, order=2,4,6,8,10 */
    	if (!sf_getfloat("phost",&phost)) phost=0;
	/* phost% points on host with zero-copy pinned memory, the rest on device */

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

	csdgather=(csd>0)?true:false;	
	_dx=1.0/dx;
	_dz=1.0/dz;
	nt_h=0.01*phost*nt; 
	nz=int((nz1+Block_Size2-1)/Block_Size2)*Block_Size2;
	nx=int((nx1+Block_Size1-1)/Block_Size1)*Block_Size1;
    	nnz=2*npml+nz; 	
	nnx=2*npml+nx;
    	dimb=dim3(Block_Size1, Block_Size2);   	dimg=dim3(nnz/Block_Size1, nnx/Block_Size2);
	dimbx1=dim3(Block_Size1, 32);		dimgx1=dim3(nnz/Block_Size1, 2);
	dimbz1=dim3(32, Block_Size2); 		dimgz1=dim3(2, nnx/Block_Size2);
	dimbx2=dim3(nz/Block_Size1,(NJ+15)/16); dimgx2=dim3(Block_Size1, 16);
	dimbz2=dim3(16, Block_Size2);		dimgz2=dim3((NJ+15)/16, nx/Block_Size2);

	if (!(dobs=(float*)malloc(ng*nt*sizeof(float)))) { sf_warning("out of memory!"); exit(1);}
	if (!(v0=(float*)malloc(nx1*nz1*sizeof(float)))) { sf_warning("out of memory!"); exit(1);}
	if (!(vel=(float*)malloc(nnz*nnx*sizeof(float)))){ sf_warning("out of memory!"); exit(1);}    	
	if (!(p=(float*)malloc(nnz*nnx*sizeof(float))))  { sf_warning("out of memory!"); exit(1);}    	
	if (!(objval=(float*)malloc(niter*sizeof(float)))){ sf_warning("out of memory!"); exit(1);}

	sf_floatread(v0, nz1*nx1, vinit);
	expand(vel, v0, npml, nnz, nnx, nz1, nx1);
    	memset(dobs, 0, ng*nt*sizeof(float));
    	memset(p, 0, nnz*nnx*sizeof(float));
	memset(objval,0,niter*sizeof(float));

    	cudaSetDevice(0);
    	cudaError_t err = cudaGetLastError();
    	if (cudaSuccess != err) 
	sf_warning("Cuda error: Failed to initialize device: %s\n", cudaGetErrorString(err));
	device_alloc();

	cuda_ricker_wavelet<<<(nt+511)/512,512>>>(d_wlt, fm, dt, nt);
	if (!(sxbeg>=0 && szbeg>=0 && sxbeg+(ns-1)*jsx<nx1 && szbeg+(ns-1)*jsz<nz1))	
	{ printf("sources exceeds the computing zone!\n"); exit(1);}
	cuda_set_sg<<<(ns+511)/512,512>>>(d_sxz, sxbeg, szbeg, jsx, jsz, ns, npml, nnz);
	distx=sxbeg-gxbeg;
	distz=szbeg-gzbeg;
	if (csdgather)	{
		if (!(gxbeg>=0 && gzbeg>=0 && gxbeg+(ng-1)*jgx<nx1 && gzbeg+(ng-1)*jgz<nz1 &&
		(sxbeg+(ns-1)*jsx)+(ng-1)*jgx-distx <nx  && (szbeg+(ns-1)*jsz)+(ng-1)*jgz-distz <nz))	
		{ sf_warning("geophones exceeds the computing zone!\n"); exit(1);}
	}
	else{
		if (!(gxbeg>=0 && gzbeg>=0 && gxbeg+(ng-1)*jgx<nx1 && gzbeg+(ng-1)*jgz<nz1))	
		{ sf_warning("geophones exceeds the computing zone!\n"); exit(1);}
	}
	cuda_set_sg<<<(ng+511)/512,512>>>(d_gxz, gxbeg, gzbeg, jgx, jgz, ng, npml, nnz);
    	cudaMemcpy(d_vv, vel, nnz*nnx*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemset(d_g0, 0, nnz*nnx*sizeof(float));	
	cudaMemset(d_g1, 0, nnz*nnx*sizeof(float));	
	cudaMemset(d_pars, 0, 4*sizeof(float));
	cudaMemset(d_bx1, 0, 2*npml*nnz*sizeof(float));
	cudaMemset(d_bz1, 0, 2*npml*nnx*sizeof(float));
	cudaMemset(d_bx2, 0, 2*npml*nnz*sizeof(float));
	cudaMemset(d_bz2, 0, 2*npml*nnx*sizeof(float));
    	cudaMemset(d_dcal, 0, ng*sizeof(float));
    	cudaMemset(d_dobs, 0, ng*nt*sizeof(float));
    	cudaMemset(d_derr, 0, ns*ng*nt*sizeof(float));
    	cudaMemset(d_g0, 0, nnz*nnx*sizeof(float));
    	cudaMemset(d_g1, 0, nnz*nnx*sizeof(float));
    	cudaMemset(d_cg, 0, nnz*nnx*sizeof(float));
	cudaMemset(d_power, 0, nnz*nnx*sizeof(float));
	cudaMemset(d_vtmp, 0, nnz*nnx*sizeof(float));
	cudaMemset(d_alpha1, 0, ng*sizeof(float));
	cudaMemset(d_alpha2, 0, ng*sizeof(float));
	cudaMemset(d_pars, 0, 4*sizeof(float));

	cudaEvent_t start,stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	for(iter=0; iter<niter; iter++)
	{
		cudaEventRecord(start);

		sf_seek(shots, 0L, SEEK_SET);
		cuda_init_abcx<<<dimgx1, dimbx1>>>(d_vv, d_bx1, d_bx2, dx, dz, dt, npml, nnz, nnx);
		cuda_init_abcz<<<dimgz1, dimbz1>>>(d_vv, d_bz1, d_bz2, dx, dz, dt, npml, nnz, nnx);
		cudaMemset(d_power, 0,	nnz*nnx*sizeof(float));
		cudaMemcpy(d_g0, d_g1, nnz*nnx*sizeof(float),cudaMemcpyDeviceToDevice);
		cudaMemset(d_g1, 0, nnz*nnx*sizeof(float));
	    	for(is=0; is<ns; is++)
	    	{ 
			sf_floatread(dobs, ng*nt, shots);
			matrix_transpose(dobs, nt, ng);
			cudaMemcpy(d_dobs, dobs, ng*nt*sizeof(float), cudaMemcpyHostToDevice);
			if (csdgather)	{
				gxbeg=sxbeg+is*jsx-distx;
				cuda_set_sg<<<(ng+511)/512, 512>>>(d_gxz, gxbeg, gzbeg, jgx, jgz, ng, npml, nnz);
			}
			cudaMemcpy(d_dobs, dobs, nt*ng*sizeof(float), cudaMemcpyHostToDevice);
			cudaMemset(h_boundary, 	0,	nt_h*2*(NJ-1)*(nx+nz)*sizeof(float));
		    	cudaMemset(d_boundary, 	0,	(nt-nt_h)*2*(NJ-1)*(nx+nz)*sizeof(float));
			wavefield_init(d_sp0, d_sp1, d_svx, d_svz, d_convpx, d_convpz, d_convvx, d_convvz);
			for(kt=0; kt<nt; kt++)
			{
			    	cuda_add_source<<<1, 1>>>(d_sp1, &d_wlt[kt], &d_sxz[is], 1, true);
				step_forward(d_vv, d_sp0, d_sp1, d_svx, d_svz, d_convvx, d_convvz, d_convpx, d_convpz, d_bx1, d_bz1, d_bx2, d_bz2);
				ptr=d_sp0; d_sp0=d_sp1; d_sp1=ptr;
	
				cuda_record<<<(ng+511)/512, 512>>>(d_sp0, d_dcal, d_gxz, ng);
				cuda_cal_residuals<<<(ng+511)/512, 512>>>(d_dcal, &d_dobs[kt*ng], &d_derr[is*nt*ng+kt*ng], ng);
					
				if(kt<nt_h)	cudaHostGetDevicePointer(&ptr, &h_boundary[kt*2*(NJ-1)*(nx+nz)], 0);
				else  		ptr=&d_boundary[(kt-nt_h)*2*(NJ-1)*(nx+nz)];
				cuda_rw_boundary_z<<<dimgz2, dimbz2>>>(ptr, 		  d_sp0, npml, nnz, nnx, NJ, false);
				cuda_rw_boundary_x<<<dimgx2, dimbx2>>>(&ptr[2*(NJ-1)*nx], d_sp0, npml, nnz, nnx, NJ, false);
			}

			// revese time loop
			ptr=d_sp0; d_sp0=d_sp1; d_sp1=ptr;// p0=p[nt]; p1=p[(nt-1)];
			wavefield_init(d_gp0, d_gp1, d_gvx, d_gvz, d_convpx, d_convpz, d_convvx, d_convvz);
			for(kt=nt-1; kt>-1; kt--)
			{
				// add residual wavefield record
				cuda_add_source<<<(ng+511)/512, 512>>>(d_gp1, &d_derr[is*nt*ng+kt*ng], d_gxz, ng, true);
				// backward time step residual wavefield
				step_forward(d_vv, d_gp0, d_gp1, d_gvx, d_gvz, d_convvx, d_convvz, d_convpx, d_convpz, d_bx1, d_bz1, d_bx2, d_bz2);
				ptr=d_gp0; d_gp0=d_gp1; d_gp1=ptr;	

				// read saved boundary
				if(kt<nt_h) 	cudaHostGetDevicePointer(&ptr, &h_boundary[kt*2*(NJ-1)*(nx+nz)], 0);
				else  		ptr=&d_boundary[(kt-nt_h)*2*(NJ-1)*(nx+nz)];
				cuda_rw_boundary_z<<<dimgz2, dimbz2>>>(ptr, 		  d_sp1, npml, nnz, nnx, NJ, true);
				cuda_rw_boundary_x<<<dimgx2, dimbx2>>>(&ptr[2*(NJ-1)*nx], d_sp1, npml, nnz, nnx, NJ, true);

				// calculate gradient
				cuda_cal_grad<<<dimg, dimb>>>(d_g1, d_sp1, d_gp1, d_power, _dz*_dz, _dx*_dx, npml, nnz, nnx);

				// backward time step source wavefield
				step_backward(d_vv, d_sp0, d_sp1, d_svx, d_svz);
				// subtract the wavelet
			    	cuda_add_source<<<1,1>>>(d_sp1, &d_wlt[kt], &d_sxz[is], 1, false);
				ptr=d_sp0; d_sp0=d_sp1; d_sp1=ptr;
			}
	    	}
		cuda_scale_grad<<<dimg, dimb>>>(d_g1, d_vv, d_power, npml, nnz, nnx);

		cuda_cal_objective<<<1, Block_Size>>>(&d_pars[0], d_derr, ns*ng*nt);
		cudaMemcpy(&obj, &d_pars[0], sizeof(float), cudaMemcpyDeviceToHost);

		cudaMemcpy(p, d_power, nnz*nnx*sizeof(float), cudaMemcpyDeviceToHost);
		window(v0, p, npml, nnz, nnx, nz1, nx1);
		sf_floatwrite(v0, nz1*nx1, illums);

		cudaMemcpy(p, d_g1, nnz*nnx*sizeof(float), cudaMemcpyDeviceToHost);
		window(v0, p, npml, nnz, nnx, nz1, nx1);
		sf_floatwrite(v0, nz1*nx1, grads);

		if (iter>0) cuda_cal_beta<<<1, Block_Size>>>(&d_pars[1], d_g0, d_g1, d_cg, nnz*nnx); 
		cudaMemcpy(&beta, &d_pars[1], sizeof(float), cudaMemcpyDeviceToHost);
		cuda_cal_conjgrad<<<dimg, dimb>>>(d_g1, d_cg, beta, npml, nnz, nnx);
		cuda_cal_epsilon<<<1, Block_Size>>>(d_vv, d_cg, &d_pars[2], nnz*nnx);
		cudaMemcpy(&epsil, &d_pars[2], sizeof(float), cudaMemcpyDeviceToHost);

		sf_seek(shots, 0L, SEEK_SET);
		cudaMemset(d_alpha1, 0, ng*sizeof(float));
		cudaMemset(d_alpha2, 0, ng*sizeof(float));
		cuda_cal_vtmp<<<dimg, dimb>>>(d_vtmp, d_vv, d_cg, epsil, npml, nnz, nnx);
		for(is=0;is<ns;is++)
		{
			sf_floatread(dobs, ng*nt, shots);
			matrix_transpose(dobs, nt, ng);
			cudaMemcpy(d_dobs, dobs, ng*nt*sizeof(float), cudaMemcpyHostToDevice);
			if (csdgather)	{
				gxbeg=sxbeg+is*jsx-distx;
				cuda_set_sg<<<(ng+511)/512, 512>>>(d_gxz, gxbeg, gzbeg, jgx, jgz, ng, npml, nnz);
			}
			wavefield_init(d_sp0, d_sp1, d_svx, d_svz, d_convpx, d_convpz, d_convvx, d_convvz);
			for(kt=0; kt<nt; kt++)
			{
			    	cuda_add_source<<<1, 1>>>(d_sp1, &d_wlt[kt], &d_sxz[is], 1, true);
				step_forward(d_vtmp, d_sp0, d_sp1, d_svx, d_svz, d_convvx, d_convvz, d_convpx, d_convpz, d_bx1, d_bz1, d_bx2, d_bz2);
				ptr=d_sp0; d_sp0=d_sp1; d_sp1=ptr;
	
				cuda_record<<<(ng+511)/512, 512>>>(d_sp0, d_dcal, d_gxz, ng);
				cuda_sum_alpha12<<<(ng+511)/512, 512>>>(d_alpha1, d_alpha2, d_dcal, &d_dobs[kt*ng], &d_derr[is*ng*nt+kt*ng], ng);
			}
		}
		cuda_cal_alpha<<<1,Block_Size>>>(&d_pars[3], d_alpha1, d_alpha2, epsil, ng);
		cudaMemcpy(&alpha, &d_pars[3], sizeof(float), cudaMemcpyDeviceToHost);

		cuda_update_vel<<<dimg,dimb>>>(d_vv, d_cg, alpha, npml, nnz, nnx);
		cudaMemcpy(p, d_vv, nnz*nnx*sizeof(float), cudaMemcpyDeviceToHost);
		window(v0, p, npml, nnz, nnx, nz1, nx1);
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

    	free(p);
	free(v0);
	free(vel);
	free(dobs);
	free(objval);
	device_free();

	return 0;
}
