/* 2D ISOTROPIC wave-equation finite-difference migration on GPU */
/*
  Copyright (C) 2012 The University of Western Australia
  
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

#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cusparse_v2.h> 
#include <cuComplex.h>
#include <cufft.h>

extern "C" {
#include <rsf.h>
}

#include "wem2d_kernels.cu"
#define a1 0.040315157
#define b1 0.873981642
#define a2 0.457289566
#define b2 0.222691983
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define MAX(x, y) (((x) > (y)) ? (x) : (y))

#define BLOCK1D 64
#define BLOCK2D 16

// ENTRY POINT
int main(int argc, char* argv[])
{
  int nw,nx,nz,nh;
  float dx,dz,ow,dw;
  int nxtap,hzero;
  bool verbose;
  float epsDSO,eps4D,epsNORM;
  
  /* I/O files */
  sf_file Fxig1 = NULL; /* Baseline XIG */
  sf_file Fxig2 = NULL; /* Monitor  XIG */
  
  sf_file Fvel1 = NULL; /* Baseline velocity file */
  sf_file Fvel2 = NULL; /* Monitor  velocity file */
  
  sf_file Fur1 = NULL; /* Baseline receiver wavefield @ nz-1 */
  sf_file Fus1 = NULL; /* Baseline  source  wavefield @ nz-1 */ 
  
  sf_file Fur2 = NULL; /* Monitor receiver wavefield @ nz-1 */
  sf_file Fus2 = NULL; /* Monitor source   wavefield @ nz-1 */ 
  
  sf_file Fgrd1 = NULL; /* Baseline Output gradient */
  sf_file Fgrd2 = NULL; /* Monitor  Output gradient */
 
  /*------------------------------------------------------------*/
  sf_init(argc,argv); /* init RSF */
  sf_axis ax,aw,az,ah; /* cube axes */
  
  /* Velocity models */
  Fvel1 = sf_input("in");   /* input velocity   */
  Fvel2 = sf_input("vel2"); /* input velocity   */
  
  /* Images */
  Fxig1 = sf_input("xig1"); sf_settype(Fxig1,SF_COMPLEX);/* input penalized image */
  Fxig2 = sf_input("xig2"); sf_settype(Fxig2,SF_COMPLEX);/* input penalized image */

  /* Wavefields */
  Fus1 = sf_input("us1"); sf_settype(Fus1,SF_COMPLEX);/* Input SWF wavefield @ nz-1 */
  Fur1 = sf_input("ur1"); sf_settype(Fur1,SF_COMPLEX);/* Input RWF wavefield @ nz-1 */
  Fus2 = sf_input("us2"); sf_settype(Fus2,SF_COMPLEX);/* Input SWF wavefield @ nz-1 */
  Fur2 = sf_input("ur2"); sf_settype(Fur2,SF_COMPLEX);/* Input RWF wavefield @ nz-1 */

  /* Gradients */
  Fgrd1 = sf_output("out"); /* Output gradient */
  Fgrd2 = sf_output("grd2"); /* Output gradient */

   /* Read in parameter */ 
  if(! sf_getint("nxtap",&nxtap)) nxtap=40;  /* TAPER size */
  if(! sf_getint("hzero",&hzero)) hzero=5;  /* Number of near offsets to zero */
  if(! sf_getbool("verbose",&verbose)) verbose=false; /* VERBOSITY flag */
  if(! sf_getfloat("epsDSO",&epsDSO)) epsDSO=1.f; /* Weighting for DSO penalty */
  if(! sf_getfloat("eps4D",&eps4D)) eps4D=0.f; /* Weighting for 4D penalty */
  if(! sf_getfloat("epsNORM",&epsNORM)) eps4D=0.00001; /* Weighting for 4D penalty */

  /* PARAMETER REPORT */
  sf_warning("nxtap: %d",nxtap);
  sf_warning("hzero: %d",hzero);
  sf_warning("epsDSO: %g",epsDSO);
  sf_warning("eps4D: %g",eps4D);
  sf_warning("epsNORM: %g",epsNORM);

  /* Get axes information */
  ax = sf_iaxa(Fur1,1);  if (verbose) sf_raxa(ax);
  aw = sf_iaxa(Fur1,2);  if (verbose) sf_raxa(aw);
  az = sf_iaxa(Fvel1,2); if (verbose) sf_raxa(az);
  ah = sf_iaxa(Fxig1,1); if (verbose) sf_raxa(ah);
  
  /* Pull out individual values */
  nx = sf_n(ax); dx = sf_d(ax);
  nw = sf_n(aw); ow = sf_o(aw); dw = sf_d(aw);
  nz = sf_n(az); dz = sf_d(az);
  nh = sf_n(ah); 
   
  /**************************************************************/
  // 
  // . . INIT GPU
  //
  /**************************************************************/
  int gpu;
  if (! sf_getint("gpu", &gpu)) gpu = 0; /* ID of the GPU to be used */
  sf_warning("using GPU #%d", gpu);
  cudaSetDevice(gpu);
  cudaDeviceSetCacheConfig(cudaFuncCachePreferL1);

  /**************************************************************/
  // 
  // . . INIT Velocity, XIG volume and GRD
  //
  /**************************************************************/
  
  /* ----------------VELOCITIES ----------------*/
  float *vel1_h=NULL; /* vel1 array (CPU) */
  float *vel2_h=NULL; /* vel2 array (CPU) */
  vel1_h = sf_floatalloc(nx*nz);
  vel2_h = sf_floatalloc(nx*nz);
  sf_floatread( vel1_h, nx*nz, Fvel1 );	
  sf_floatread( vel2_h, nx*nz, Fvel2 );	
  for (int ixz=0; ixz < nx*nz; ixz++) {
    vel1_h[ixz] = 1.f / vel1_h[ixz] ;
  }
  for (int ixz=0; ixz < nx*nz; ixz++) {
    vel2_h[ixz] = 1.f / vel2_h[ixz] ;
  }
  
  float *vel1_d; /* vel array (GPU) */ 
  float *vel2_d; /* vel array (GPU) */ 
  cudaMalloc((void **)&vel1_d, nx*nz*sizeof(float) );
  cudaMalloc((void **)&vel2_d, nx*nz*sizeof(float) );
  cudaMemcpy(vel1_d, vel1_h, nx*nz*sizeof(float), cudaMemcpyHostToDevice );
  cudaMemcpy(vel2_d, vel2_h, nx*nz*sizeof(float), cudaMemcpyHostToDevice );
  
  /*---------------- XIGS ----------------*/
  sf_complex *xig1_h=NULL; /* XIG (CPU) 3D */
  sf_complex *xig2_h=NULL; /* XIG (CPU) 3D */
  xig1_h = sf_complexalloc( nh * nx * nz);
  xig2_h = sf_complexalloc( nh * nx * nz);
  sf_complexread( xig1_h, nh*nx*nz, Fxig1);	
  sf_complexread( xig2_h, nh*nx*nz, Fxig2);	

  sf_complex *xigslice1_h=NULL; /* XIG (CPU) 2D slice*/
  sf_complex *xigslice2_h=NULL; /* XIG (CPU) 2D slice*/
  xigslice1_h = sf_complexalloc( nh * nx );
  xigslice2_h = sf_complexalloc( nh * nx );

  cuComplex *xig1_d; /* XIG array (GPU) */
  cuComplex *xig2_d; /* XIG array (GPU) */
  cudaMalloc((void **)&xig1_d, nx*nh*sizeof(cuComplex) );
  cudaMalloc((void **)&xig2_d, nx*nh*sizeof(cuComplex) );
  cudaMemset(xig1_d, 0, nx*nw*sizeof(cuComplex));
  cudaMemset(xig2_d, 0, nx*nw*sizeof(cuComplex));

  /* ---------------- GRADIENTS ----------------*/
  float *grd1_h=NULL; /* Gradient (2D) CPU */
  float *grd2_h=NULL; /* Gradient (2D) CPU */
  grd1_h = sf_floatalloc( nx * nz);
  grd2_h = sf_floatalloc( nx * nz);

  float *grdslice1_h=NULL; /* Gradient slice (1D) CPU */
  float *grdslice2_h=NULL; /* Gradient slice (1D) CPU */
  grdslice1_h = sf_floatalloc(nx);
  grdslice2_h = sf_floatalloc(nx);

  float *grd1_d; /* Gradient (2D) GPU */
  float *grd2_d; /* Gradient (2D) GPU */
  cudaMalloc((void **)&grd1_d,nx*sizeof(float) );
  cudaMalloc((void **)&grd2_d,nx*sizeof(float) );
  cudaMemset(grd1_d,0,nx*sizeof(float));
  cudaMemset(grd2_d,0,nx*sizeof(float));

  /**************************************************************/
  // 
  // . . Read in Wavefields from iz=nz+1
  //
  /**************************************************************/

  /* ---------------- BASELINE ----------------*/
  sf_complex *us1_h=NULL; /* SWF on CPU (2D) */
  us1_h = sf_complexalloc( nx * nw );
  sf_complexread(us1_h, nx*nw, Fus1);

  cuComplex *us1_d; /* SWF on GPU (2D) */
  cudaMalloc((void **)&us1_d,nx*nw*sizeof(cuComplex) );
  cudaMemcpy( us1_d, us1_h,  nx*nw*sizeof(cuComplex), cudaMemcpyHostToDevice );
  
  sf_complex *ur1_h=NULL; /* RWF on CPU (2D) */
  ur1_h = sf_complexalloc( nx * nw );
  sf_complexread(ur1_h, nx*nw, Fur1);

  cuComplex *ur1_d; /* RWF on GPU (2D) */
  cudaMalloc((void **)&ur1_d,nx*nw*sizeof(cuComplex) );
  cudaMemcpy( ur1_d, ur1_h,  nx*nw*sizeof(cuComplex), cudaMemcpyHostToDevice );
  
  /* ---------------- MONITOR ----------------*/
  sf_complex *us2_h=NULL; /* SWF on CPU (2D) */
  us2_h = sf_complexalloc( nx * nw );
  sf_complexread(us2_h, nx*nw, Fus2);

  cuComplex *us2_d; /* SWF on GPU (2D) */
  cudaMalloc((void **)&us2_d,nx*nw*sizeof(cuComplex) );
  cudaMemcpy( us2_d, us2_h,  nx*nw*sizeof(cuComplex), cudaMemcpyHostToDevice );
  
  sf_complex *ur2_h=NULL; /* RWF on CPU (2D) */
  ur2_h = sf_complexalloc( nx * nw );
  sf_complexread(ur2_h, nx*nw, Fur2);

  cuComplex *ur2_d; /* RWF on GPU (2D) */
  cudaMalloc((void **)&ur2_d,nx*nw*sizeof(cuComplex) );
  cudaMemcpy( ur2_d, ur2_h,  nx*nw*sizeof(cuComplex), cudaMemcpyHostToDevice );
  
  /* Initialize Adjoint wavefield volumes on GPUs and set to zero */
  /* ---------------- BASELINE ----------------*/
  cuComplex *as1_d;
  cudaMalloc((void **)&as1_d, nx*nw*sizeof(cuComplex)); 
  cudaMemset(as1_d, 0, nx*nw*sizeof(cuComplex));

  cuComplex *ar1_d;
  cudaMalloc((void **)&ar1_d, nx*nw*sizeof(cuComplex)); 
  cudaMemset(ar1_d, 0, nx*nw*sizeof(cuComplex));

  /* ---------------- MONITOR ----------------*/
  cuComplex *as2_d;
  cudaMalloc((void **)&as2_d, nx*nw*sizeof(cuComplex)); 
  cudaMemset(as2_d, 0, nx*nw*sizeof(cuComplex));

  cuComplex *ar2_d;
  cudaMalloc((void **)&ar2_d, nx*nw*sizeof(cuComplex)); 
  cudaMemset(ar2_d, 0, nx*nw*sizeof(cuComplex));

  /**************************************************************/
  /* Initialize wavefield taper and copy to GPU */
  sf_complex *tap_h=NULL; /* taper on CPU */
  cuComplex *tap_d; /* taper on GPU */
  tap_h = sf_complexalloc(nx);
  int j1;
  for (int ix=0; ix<nx; ix++) {
  	tap_h[ix]= sf_cmplx(1.f,0.f);
  }
  for (int ix=0; ix < nxtap; ix++) {
  	j1 = abs(nxtap-ix-1.f);
    tap_h[ix] = sf_cmplx( cos ( SF_PI/2.f*j1/(nxtap-1.f) ), 0.f);
  }
  
  for (int ix=0; ix < nxtap; ix++) {
  	j1 = abs(nxtap-ix-1.f);
    tap_h[nx-ix-1] = sf_cmplx( cos (SF_PI/2.f*j1/(nxtap-1.f) ),0.f);
  }
    
  tap_h[0]=sf_cmplx(0.f,0.f); 
  tap_h[nx-1]=sf_cmplx(0.f,0.f);
  
  cudaMalloc((void **)&tap_d, nx*sizeof(cuComplex) );
  cudaMemcpy( tap_d, tap_h,   nx*sizeof(cuComplex), cudaMemcpyHostToDevice );

  /**************************************************************/
  /* Set up vmin array */
  float *vmin1_h=NULL;
  float *vmin2_h=NULL;
  vmin1_h = sf_floatalloc(nz);
  vmin2_h = sf_floatalloc(nz);
  float m;
  for (int  iz=0; iz<nz; iz++) {
	m = SF_HUGE;  
	for (int ix=0; ix<nx; ix++) {
		if (vel1_h[iz*nx+ix] < m) m = vel1_h[iz*nx+ix];
	}
	vmin1_h[iz]=m;
  }	
  for (int  iz=0; iz<nz; iz++) {
	m = SF_HUGE;  
	for (int ix=0; ix<nx; ix++) {
		if (vel2_h[iz*nx+ix] < m) m = vel2_h[iz*nx+ix];
	}
	vmin2_h[iz]=m;
  }	
  /**************************************************************/
  /* Set up Fourier domain wavenumbers */	
  float *kxmap_h=NULL; float *kxmap_d;
  kxmap_h = sf_floatalloc(nx); 
  int mm=nx/2;
  float kscale = 2.f*SF_PI/( (float)nx * dx);
  for (int ii=0; ii<mm+1; ii++) {
  	kxmap_h[ii]=kscale*((float)ii);
  }	

  for (int ii=mm; ii<nx; ii++) {
  	kxmap_h[ii]=kscale*((float)ii-(float)nx);
  }
  cudaMalloc((void **)&kxmap_d, nx*sizeof(float) );
  cudaMemcpy( kxmap_d, kxmap_h, nx*sizeof(float), cudaMemcpyHostToDevice );
 
 
  /**************************************************************/
  /* Set aabb array FD coeffs including dt and dx */
  float *aabb_h=NULL;
  aabb_h = sf_floatalloc(4);
  aabb_h[0] = a1 * dz / (2.f*dx*dx);
  aabb_h[1] = b1 / (dx * dx);
  aabb_h[2] = a2 * dz / (2.f*dx*dx);
  aabb_h[3] = b2 / (dx * dx);
  
  /**************************************************************/
  /* Set up Array calls */	
  dim3 dimBlock(BLOCK2D,BLOCK2D,1); 
  dim3 dimGrid( ceil(nx/(float)BLOCK2D), ceil(nw/(float)BLOCK2D) ,1); 

  int imgcall=BLOCK1D; float imgcallf=(float)imgcall;
  dim3 dimBlock2( imgcall, 1);
  dim3 dimGrid2( ceil(nx/imgcallf), 1 );

  /* Define causality */
  float  caus= 1.f;
  float acaus=-1.f;

  /**************************************************************/
  /* Set up arrays required by FD code */	
  cuComplex *rax_d, *rbx_d, *rcx_d, *v_d;
  
  /* Allocate Memory on GPUs */
  cudaMalloc((void **)&rax_d, nx*nw*sizeof(cuComplex)); 
  cudaMalloc((void **)&rbx_d, nx*nw*sizeof(cuComplex)); 
  cudaMalloc((void **)&rcx_d, nx*nw*sizeof(cuComplex)); 
  cudaMalloc((void **)&v_d,   nx*nw*sizeof(cuComplex)); 

  /* Set equal to zero */
  cudaMemset(rax_d, 0, nx*nw*sizeof(cuComplex));
  cudaMemset(rbx_d, 0, nx*nw*sizeof(cuComplex));
  cudaMemset(rcx_d, 0, nx*nw*sizeof(cuComplex));
  cudaMemset(v_d,   0, nx*nw*sizeof(cuComplex));

  /**************************************************************/
  /* FFT set up Create a 1D FFT plan. */
  cufftHandle plan;
  cufftPlan1d(&plan, nx, CUFFT_C2C, nw);

  /**************************************************************/
  /* Handle tridiagonal solver */
  cusparseHandle_t cusparseHandle = 0;
//  cusparseStatus_t cusparseStatus;
//  cusparseStatus = cusparseCreate(&cusparseHandle);

  /**************************************************************/
  /**************************************************************/
  /**************************************************************/
  /*                        MAIN LOOP                           */	
  /**************************************************************/
  /**************************************************************/
  /**************************************************************/

  /**************************************************************/
  // 
  // . . PASS ON XIG TO DEVICE @ iz=nz-1
  //
  /**************************************************************/	
  
  /* BASELINE */
  for (int ixh=0; ixh < nh*nx; ixh++) {xigslice1_h[ixh] = xig1_h[(nz-1)*nx*nh+ixh];}
  cudaMemcpy(xig1_d, xigslice1_h, nx*nh*sizeof(cuComplex), cudaMemcpyHostToDevice );
  
  /* MONITOR */
  for (int ixh=0; ixh < nh*nx; ixh++) {xigslice2_h[ixh] = xig2_h[(nz-1)*nx*nh+ixh];}
  cudaMemcpy(xig2_d, xigslice2_h, nx*nh*sizeof(cuComplex), cudaMemcpyHostToDevice );

  /**************************************************************/
  // 
  // . . COMPUTE ADJOINT WAVEFIELDS (COUPLED)
  //
  /**************************************************************/	
  make_adj_wflds_coupled<<<dimGrid,dimBlock>>> (xig1_d,xig2_d,us1_d,ur1_d,as1_d,ar1_d,us2_d,ur2_d,as2_d,ar2_d,epsDSO,hzero,eps4D,epsNORM,dx,nh,nx,nw);

  /**************************************************************/
  // 
  // . . COMPUTE GRADIENT
  //
  /**************************************************************/	
  compute_gradient<<<dimGrid2,dimBlock2>>>(grd1_d,us1_d,ur1_d,as1_d,ar1_d,nx,ow,dw,nw);	 /* BASELINE */	
  compute_gradient<<<dimGrid2,dimBlock2>>>(grd2_d,us2_d,ur2_d,as2_d,ar2_d,nx,ow,dw,nw);	 /* MONITOR  */		

  /**************************************************************/
  // 
  // . . COPY BACK GRADIENT TO HOST
  //
  /**************************************************************/	
  /* BASELINE */	
  cudaMemcpy( grdslice1_h, grd1_d, nx*sizeof(float), cudaMemcpyDeviceToHost );
  for (int ix=0; ix<nx; ix++) { grd1_h[(nz-1)*nx+ix]=grdslice1_h[ix];}

  /* MONITOR  */	
  cudaMemcpy( grdslice2_h, grd2_d, nx*sizeof(float), cudaMemcpyDeviceToHost );
  for (int ix=0; ix<nx; ix++) { grd2_h[(nz-1)*nx+ix]=grdslice2_h[ix];}
  
  /**************************************************************/
  // 
  // . . STEP FROM iz=NZ-2 to 0
  //
  /**************************************************************/	  
  for (int iz = nz-2; iz > -1; iz--) {  /* Depth Loop */ 
    	if (verbose) fprintf(stderr,"\b\b\b\b\b\b\b\b\b\b\b\b%d of %d",nz-iz,nz); 

		/* EXTRACT image slice */
  		for (int ixh=0; ixh < nh*nx; ixh++) { xigslice1_h[ixh] = xig1_h[iz*nx*nh+ixh];}
  		for (int ixh=0; ixh < nh*nx; ixh++) { xigslice2_h[ixh] = xig2_h[iz*nx*nh+ixh];}
		
		/* Copy from host to device */
		cudaMemcpy(xig1_d,xigslice1_h,nx*nh*sizeof(cuComplex), cudaMemcpyHostToDevice );
  		cudaMemcpy(xig2_d,xigslice2_h,nx*nh*sizeof(cuComplex), cudaMemcpyHostToDevice );
  	
   		/**************************************************************/
  		// 
  		// . . WAVEFIELD TAPER
  		//
  		/**************************************************************/
  		
  		/* BASELINE */
  		wfld_taper<<<dimGrid, dimBlock >>>(us1_d,tap_d,nx,nw);
		wfld_taper<<<dimGrid, dimBlock >>>(ur1_d,tap_d,nx,nw);
  		wfld_taper<<<dimGrid, dimBlock >>>(as1_d,tap_d,nx,nw);
		wfld_taper<<<dimGrid, dimBlock >>>(ar1_d,tap_d,nx,nw);
		
		/* MONITOR */
  		wfld_taper<<<dimGrid, dimBlock >>>(us2_d,tap_d,nx,nw);
		wfld_taper<<<dimGrid, dimBlock >>>(ur2_d,tap_d,nx,nw);
  		wfld_taper<<<dimGrid, dimBlock >>>(as2_d,tap_d,nx,nw);
		wfld_taper<<<dimGrid, dimBlock >>>(ar2_d,tap_d,nx,nw);

   		/**************************************************************/
  		// 
  		// . . SPLIT-STEP PROPAGATION
  		//
  		/**************************************************************/
  		
  		/* BASELINE */
  		prop_SSF<<<dimGrid,dimBlock>>>(us1_d,vel1_d,dw,ow,acaus,iz,dz,nx,nw);
		prop_SSF<<<dimGrid,dimBlock>>>(ur1_d,vel1_d,dw,ow, caus,iz,dz,nx,nw);
  		prop_SSF<<<dimGrid,dimBlock>>>(as1_d,vel1_d,dw,ow,acaus,iz,dz,nx,nw);
		prop_SSF<<<dimGrid,dimBlock>>>(ar1_d,vel1_d,dw,ow, caus,iz,dz,nx,nw);
		
		/* MONITOR */
  		prop_SSF<<<dimGrid,dimBlock>>>(us2_d,vel2_d,dw,ow,acaus,iz,dz,nx,nw);
		prop_SSF<<<dimGrid,dimBlock>>>(ur2_d,vel2_d,dw,ow, caus,iz,dz,nx,nw);
  		prop_SSF<<<dimGrid,dimBlock>>>(as2_d,vel2_d,dw,ow,acaus,iz,dz,nx,nw);
		prop_SSF<<<dimGrid,dimBlock>>>(ar2_d,vel2_d,dw,ow, caus,iz,dz,nx,nw);

  		/**************************************************************/
  		// 
  		// . . STEP WAVEFIELDS using tridiagonal and CUSPARSE
  		//
  		/**************************************************************/
  		
  		/* BASELINE */
		setup_FD<<<dimGrid,dimBlock>>>(us1_d,v_d,rax_d,rbx_d,rcx_d,vel1_d,dw,ow,-acaus*aabb_h[0],aabb_h[1],nx,nw,iz);
		cusparseCgtsvStridedBatch(cusparseHandle, nx, rcx_d, rax_d, rbx_d, v_d, nw, nx);
		copy_wfld<<<dimGrid,dimBlock>>>(v_d,us1_d,nx,nw); /* US1 */
		setup_FD<<<dimGrid,dimBlock>>>(us1_d,v_d,rax_d,rbx_d,rcx_d,vel1_d,dw,ow,-acaus*aabb_h[2],aabb_h[3],nx,nw,iz);
		cusparseCgtsvStridedBatch(cusparseHandle, nx, rcx_d, rax_d, rbx_d, v_d, nw, nx);
		copy_wfld<<<dimGrid,dimBlock>>>(v_d,us1_d,nx,nw);

		setup_FD<<<dimGrid,dimBlock>>>(ur1_d,v_d,rax_d,rbx_d,rcx_d,vel1_d,dw,ow, -caus*aabb_h[0],aabb_h[1],nx,nw,iz);
		cusparseCgtsvStridedBatch(cusparseHandle, nx, rcx_d, rax_d, rbx_d, v_d, nw, nx);
		copy_wfld<<<dimGrid,dimBlock>>>(v_d,ur1_d,nx,nw); /* UR1 */
		setup_FD<<<dimGrid,dimBlock>>>(ur1_d,v_d,rax_d,rbx_d,rcx_d,vel1_d,dw,ow, -caus*aabb_h[2],aabb_h[3],nx,nw,iz);
		cusparseCgtsvStridedBatch(cusparseHandle, nx, rcx_d, rax_d, rbx_d, v_d, nw, nx);
		copy_wfld<<<dimGrid,dimBlock>>>(v_d,ur1_d,nx,nw);
 
		setup_FD<<<dimGrid,dimBlock>>>(as1_d,v_d,rax_d,rbx_d,rcx_d,vel1_d,dw,ow,-acaus*aabb_h[0],aabb_h[1],nx,nw,iz);
		cusparseCgtsvStridedBatch(cusparseHandle, nx, rcx_d, rax_d, rbx_d, v_d, nw, nx);
		copy_wfld<<<dimGrid,dimBlock>>>(v_d,as1_d,nx,nw); /* AS1 */
		setup_FD<<<dimGrid,dimBlock>>>(as1_d,v_d,rax_d,rbx_d,rcx_d,vel1_d,dw,ow,-acaus*aabb_h[2],aabb_h[3],nx,nw,iz);
		cusparseCgtsvStridedBatch(cusparseHandle, nx, rcx_d, rax_d, rbx_d, v_d, nw, nx);
		copy_wfld<<<dimGrid,dimBlock>>>(v_d,as1_d,nx,nw);

		setup_FD<<<dimGrid,dimBlock>>>(ar1_d,v_d,rax_d,rbx_d,rcx_d,vel1_d,dw,ow, -caus*aabb_h[0],aabb_h[1],nx,nw,iz);
		cusparseCgtsvStridedBatch(cusparseHandle, nx, rcx_d, rax_d, rbx_d, v_d, nw, nx);
		copy_wfld<<<dimGrid,dimBlock>>>(v_d,ar1_d,nx,nw); /* AS1 */
		setup_FD<<<dimGrid,dimBlock>>>(ar1_d,v_d,rax_d,rbx_d,rcx_d,vel1_d,dw,ow, -caus*aabb_h[2],aabb_h[3],nx,nw,iz);
		cusparseCgtsvStridedBatch(cusparseHandle, nx, rcx_d, rax_d, rbx_d, v_d, nw, nx);
		copy_wfld<<<dimGrid,dimBlock>>>(v_d,ar1_d,nx,nw);
  
  	 	/* MONITOR */
		setup_FD<<<dimGrid,dimBlock>>>(us2_d,v_d,rax_d,rbx_d,rcx_d,vel2_d,dw,ow,-acaus*aabb_h[0],aabb_h[1],nx,nw,iz);
		cusparseCgtsvStridedBatch(cusparseHandle, nx, rcx_d, rax_d, rbx_d, v_d, nw, nx);
		copy_wfld<<<dimGrid,dimBlock>>>(v_d,us2_d,nx,nw); /* US2 */
		setup_FD<<<dimGrid,dimBlock>>>(us2_d,v_d,rax_d,rbx_d,rcx_d,vel2_d,dw,ow,-acaus*aabb_h[2],aabb_h[3],nx,nw,iz);
		cusparseCgtsvStridedBatch(cusparseHandle, nx, rcx_d, rax_d, rbx_d, v_d, nw, nx);
		copy_wfld<<<dimGrid,dimBlock>>>(v_d,us2_d,nx,nw);

		setup_FD<<<dimGrid,dimBlock>>>(ur2_d,v_d,rax_d,rbx_d,rcx_d,vel2_d,dw,ow, -caus*aabb_h[0],aabb_h[1],nx,nw,iz);
		cusparseCgtsvStridedBatch(cusparseHandle, nx, rcx_d, rax_d, rbx_d, v_d, nw, nx);
		copy_wfld<<<dimGrid,dimBlock>>>(v_d,ur2_d,nx,nw); /* UR2 */
		setup_FD<<<dimGrid,dimBlock>>>(ur2_d,v_d,rax_d,rbx_d,rcx_d,vel2_d,dw,ow, -caus*aabb_h[2],aabb_h[3],nx,nw,iz);
		cusparseCgtsvStridedBatch(cusparseHandle, nx, rcx_d, rax_d, rbx_d, v_d, nw, nx);
		copy_wfld<<<dimGrid,dimBlock>>>(v_d,ur2_d,nx,nw);
 
		setup_FD<<<dimGrid,dimBlock>>>(as2_d,v_d,rax_d,rbx_d,rcx_d,vel2_d,dw,ow,-acaus*aabb_h[0],aabb_h[1],nx,nw,iz);
		cusparseCgtsvStridedBatch(cusparseHandle, nx, rcx_d, rax_d, rbx_d, v_d, nw, nx);
		copy_wfld<<<dimGrid,dimBlock>>>(v_d,as2_d,nx,nw); /* AS2 */
		setup_FD<<<dimGrid,dimBlock>>>(as2_d,v_d,rax_d,rbx_d,rcx_d,vel2_d,dw,ow,-acaus*aabb_h[2],aabb_h[3],nx,nw,iz);
		cusparseCgtsvStridedBatch(cusparseHandle, nx, rcx_d, rax_d, rbx_d, v_d, nw, nx);
		copy_wfld<<<dimGrid,dimBlock>>>(v_d,as2_d,nx,nw);

		setup_FD<<<dimGrid,dimBlock>>>(ar2_d,v_d,rax_d,rbx_d,rcx_d,vel2_d,dw,ow, -caus*aabb_h[0],aabb_h[1],nx,nw,iz);
		cusparseCgtsvStridedBatch(cusparseHandle, nx, rcx_d, rax_d, rbx_d, v_d, nw, nx);
		copy_wfld<<<dimGrid,dimBlock>>>(v_d,ar2_d,nx,nw); /* AS2 */
		setup_FD<<<dimGrid,dimBlock>>>(ar2_d,v_d,rax_d,rbx_d,rcx_d,vel2_d,dw,ow, -caus*aabb_h[2],aabb_h[3],nx,nw,iz);
		cusparseCgtsvStridedBatch(cusparseHandle, nx, rcx_d, rax_d, rbx_d, v_d, nw, nx);
		copy_wfld<<<dimGrid,dimBlock>>>(v_d,ar2_d,nx,nw);  
  
   		/**************************************************************/
  		// 
  		// . . FORWARD FOURIER TRANSFORM
  		//
  		/**************************************************************/
  		
  		/* BASELINE */
 		cufftExecC2C(plan,us1_d,us1_d,CUFFT_FORWARD);
		cufftExecC2C(plan,ur1_d,ur1_d,CUFFT_FORWARD);
		cufftExecC2C(plan,as1_d,as1_d,CUFFT_FORWARD);
		cufftExecC2C(plan,ar1_d,ar1_d,CUFFT_FORWARD);
		
		/* MONITOR */
 		cufftExecC2C(plan,us2_d,us2_d,CUFFT_FORWARD);
		cufftExecC2C(plan,ur2_d,ur2_d,CUFFT_FORWARD);
		cufftExecC2C(plan,as2_d,as2_d,CUFFT_FORWARD);
		cufftExecC2C(plan,ar2_d,ar2_d,CUFFT_FORWARD);
		
   		/**************************************************************/
  		// 
  		// . . PHASE CORRECTIONS
  		//
  		/**************************************************************/
  		
  		/* BASELINE */
 		phase_correction<<<dimGrid,dimBlock>>>(us1_d,vmin1_h[iz],dw,ow,kxmap_d,iz,-acaus,nx,nw,dz);
		phase_correction<<<dimGrid,dimBlock>>>(ur1_d,vmin1_h[iz],dw,ow,kxmap_d,iz,- caus,nx,nw,dz);
 		phase_correction<<<dimGrid,dimBlock>>>(as1_d,vmin1_h[iz],dw,ow,kxmap_d,iz,-acaus,nx,nw,dz);
		phase_correction<<<dimGrid,dimBlock>>>(ar1_d,vmin1_h[iz],dw,ow,kxmap_d,iz,- caus,nx,nw,dz);
		
		/* MONITOR */
 		phase_correction<<<dimGrid,dimBlock>>>(us2_d,vmin2_h[iz],dw,ow,kxmap_d,iz,-acaus,nx,nw,dz);
		phase_correction<<<dimGrid,dimBlock>>>(ur2_d,vmin2_h[iz],dw,ow,kxmap_d,iz,- caus,nx,nw,dz);
 		phase_correction<<<dimGrid,dimBlock>>>(as2_d,vmin2_h[iz],dw,ow,kxmap_d,iz,-acaus,nx,nw,dz);
		phase_correction<<<dimGrid,dimBlock>>>(ar2_d,vmin2_h[iz],dw,ow,kxmap_d,iz,- caus,nx,nw,dz);
	
   		/**************************************************************/
  		// 
  		// . . INVERSE FOURIER TRANSFORM
  		//
  		/**************************************************************/
  		
  		/* BASELINE */
  		cufftExecC2C(plan,us1_d,us1_d,CUFFT_INVERSE);
		cufftExecC2C(plan,ur1_d,ur1_d,CUFFT_INVERSE);
  		cufftExecC2C(plan,as1_d,as1_d,CUFFT_INVERSE);
		cufftExecC2C(plan,ar1_d,ar1_d,CUFFT_INVERSE);
		
		/* MONITOR */
  		cufftExecC2C(plan,us2_d,us2_d,CUFFT_INVERSE);
		cufftExecC2C(plan,ur2_d,ur2_d,CUFFT_INVERSE);
  		cufftExecC2C(plan,as2_d,as2_d,CUFFT_INVERSE);
		cufftExecC2C(plan,ar2_d,ar2_d,CUFFT_INVERSE);
		
   		/**************************************************************/
  		// 
  		// . . RESCALE WAVEFIELDS
  		//
  		/**************************************************************/
  		
  		/* BASELINE */
  		rescale_wfld<<< dimGrid,dimBlock >>>(us1_d,nx,nw);
		rescale_wfld<<< dimGrid,dimBlock >>>(ur1_d,nx,nw);
  		rescale_wfld<<< dimGrid,dimBlock >>>(as1_d,nx,nw);
		rescale_wfld<<< dimGrid,dimBlock >>>(ar1_d,nx,nw);
		
  		/* MONITOR */
  		rescale_wfld<<< dimGrid,dimBlock >>>(us2_d,nx,nw);
		rescale_wfld<<< dimGrid,dimBlock >>>(ur2_d,nx,nw);
  		rescale_wfld<<< dimGrid,dimBlock >>>(as2_d,nx,nw);
		rescale_wfld<<< dimGrid,dimBlock >>>(ar2_d,nx,nw);
	
   		/**************************************************************/
  		// 
  		// . . WAVEFIELD TAPER
  		//
  		/**************************************************************/
  		
  		/* BASELINE */
  		wfld_taper<<<dimGrid,dimBlock>>>(us1_d,tap_d,nx,nw);
		wfld_taper<<<dimGrid,dimBlock>>>(ur1_d,tap_d,nx,nw);
  		wfld_taper<<<dimGrid,dimBlock>>>(as1_d,tap_d,nx,nw);
		wfld_taper<<<dimGrid,dimBlock>>>(ar1_d,tap_d,nx,nw);
	
		/* MONITOR */
  		wfld_taper<<<dimGrid,dimBlock>>>(us2_d,tap_d,nx,nw);
		wfld_taper<<<dimGrid,dimBlock>>>(ur2_d,tap_d,nx,nw);
  		wfld_taper<<<dimGrid,dimBlock>>>(as2_d,tap_d,nx,nw);
		wfld_taper<<<dimGrid,dimBlock>>>(ar2_d,tap_d,nx,nw);
		
   		/**************************************************************/
  		// 
  		// . . COMPUTE ADJOINT WAVEFIELDS
  		//
  		/**************************************************************/	
		
  		/* Make the adjoint wavefields */
 		make_adj_wflds_coupled<<<dimGrid,dimBlock>>> (xig1_d,xig2_d,us1_d,ur1_d,as1_d,ar1_d,us2_d,ur2_d,as2_d,ar2_d,epsDSO,hzero,eps4D,epsNORM,dx,nh,nx,nw);
	
   		/**************************************************************/
  		// 
  		// . . COMPUTE GRADIENT 
  		//
  		/**************************************************************/	
  		compute_gradient<<<dimGrid2,dimBlock2>>>(grd1_d,us1_d,ur1_d,as1_d,ar1_d,nx,ow,dw,nw); /* BASELINE */		
  		compute_gradient<<<dimGrid2,dimBlock2>>>(grd2_d,us2_d,ur2_d,as2_d,ar2_d,nx,ow,dw,nw); /* MONITOR  */

		/* Copy back XIG slice and write output (include sum over frequency) */
		/* BASELINE */	
  		cudaMemcpy( grdslice1_h, grd1_d, nx*sizeof(float), cudaMemcpyDeviceToHost );
  		for (int ix=0; ix<nx; ix++) { grd1_h[iz*nx+ix]=grdslice1_h[ix];} 
  		/* MONITOR */
   		cudaMemcpy( grdslice2_h, grd2_d, nx*sizeof(float), cudaMemcpyDeviceToHost );
  		for (int ix=0; ix<nx; ix++) { grd2_h[iz*nx+ix]=grdslice2_h[ix];} 		
  		
    } /* END DEPTH LOOP */

  /* Write out gradient */
  if (verbose) fprintf(stderr,"\n"); 
  sf_warning("WRITING OUTPUT GRADIENTS");
  
  sf_floatwrite(grd1_h, nx*nz, Fgrd1);
  sf_floatwrite(grd2_h, nx*nz, Fgrd2);

  /**************************************************************/
  // 
  // . . COMPUTE GRADIENT 
  //
  /**************************************************************/
  cufftDestroy(plan);
  cudaFree(vel1_d); cudaFree(xig1_d); cudaFree(grd1_d);   
  cudaFree(us1_d);  cudaFree(ur1_d);
  cudaFree(as1_d);  cudaFree(ar1_d);

  cudaFree(vel2_d); cudaFree(xig2_d); cudaFree(grd2_d);   
  cudaFree(us2_d);  cudaFree(ur2_d);
  cudaFree(as2_d);  cudaFree(ar2_d);
    
  cudaFree(rax_d); cudaFree(rbx_d); cudaFree(rcx_d);
  cudaFree(v_d); cudaFree(tap_d);   cudaFree(kxmap_d);
  

  exit(0);
}
