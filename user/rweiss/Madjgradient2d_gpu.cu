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
  int nxtap;
  bool verbose;
#if (CUDART_VERSION >= 10000)
  size_t pbuffersize;
#endif
  
  /* I/O files */
  sf_file Fxig = NULL; /* Input XIG */
  sf_file Fvel = NULL; /* velocity file */
  sf_file Frwf = NULL; /* receiver wavefield @ nz-1 */
  sf_file Fswf = NULL; /* source   wavefield @ nz-1 */ 
  sf_file Fgrd = NULL; /* Output gradient */
 
  /*------------------------------------------------------------*/
  sf_init(argc,argv); /* init RSF */
  sf_axis ax,aw,az,ah; /* cube axes */
  
  /* Set up file objects */
  Fvel = sf_input("in");  /* input velocity   */
  Fxig = sf_input("xig"); sf_settype(Fxig,SF_COMPLEX);/* input penalized image */
  Fswf = sf_input("swf"); sf_settype(Fswf,SF_COMPLEX);/* Input SWF wavefield @ nz-1 */
  Frwf = sf_input("rwf"); sf_settype(Frwf,SF_COMPLEX);/* Input RWF wavefield @ nz-1 */
  Fgrd = sf_output("out"); /* Output gradient */

   /* Read in parameter */ 
  if(! sf_getint("nxtap",&nxtap)) nxtap=40;  /* TAPER size */
  if(! sf_getbool("verbose",&verbose)) verbose=false; /* VERBOSITY flag */

  /* Get axes information */
  ax = sf_iaxa(Frwf,1); if (verbose) sf_raxa(ax);
  az = sf_iaxa(Fvel,2); if (verbose) sf_raxa(az);
  aw = sf_iaxa(Frwf,2); if (verbose) sf_raxa(aw);
  ah = sf_iaxa(Fxig,1); if (verbose) sf_raxa(ah);
  
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
  float *vel_h=NULL; /* vel array (CPU) */
  vel_h = sf_floatalloc(nx*nz);
  sf_floatread( vel_h, nx*nz, Fvel );	
  for (int ixz=0; ixz < nx*nz; ixz++) {
    vel_h[ixz] = 1.f / vel_h[ixz] ;
  }

  float *vel_d; /* vel array (GPU) */ 
  cudaMalloc((void **)&vel_d, nx*nz*sizeof(float) );
  cudaMemcpy(vel_d, vel_h, nx*nz*sizeof(float), cudaMemcpyHostToDevice );
  
  sf_complex *xig_h=NULL; /* XIG (CPU) 3D */
  xig_h = sf_complexalloc( nh * nx * nz);
  sf_complexread( xig_h, nh*nx*nz, Fxig );	

  sf_complex *xigslice_h=NULL; /* XIG (CPU) 2D slice*/
  xigslice_h = sf_complexalloc( nh * nx );

  cuComplex *xig_d; /* XIG array (GPU) */
  cudaMalloc((void **)&xig_d, nx*nh*sizeof(cuComplex) );
  cudaMemset(xig_d, 0, nx*nw*sizeof(cuComplex));

  float *grd_h=NULL; /* Gradient (2D) CPU */
  grd_h = sf_floatalloc( nx * nz);

  float *grdslice_h=NULL; /* Gradient slice (1D) CPU */
  grdslice_h = sf_floatalloc(nx);

  float *grd_d; /* Gradient (2D) GPU */
  cudaMalloc((void **)&grd_d, nx*sizeof(float) );
  cudaMemset(grd_d, 0, nx*sizeof(cuComplex));

  /**************************************************************/
  // 
  // . . Read in Wavefields from iz=nz+1
  //
  /**************************************************************/
  sf_complex *swf_h=NULL; /* SWF on CPU (2D) */
  swf_h = sf_complexalloc( nx * nw );
  sf_complexread(swf_h, nx*nw, Fswf);

  cuComplex *swf_d; /* SWF on GPU (2D) */
  cudaMalloc((void **)&swf_d, nx*nw*sizeof(cuComplex) );
  cudaMemcpy( swf_d, swf_h, nx*nw*sizeof(cuComplex), cudaMemcpyHostToDevice );
  
  sf_complex *rwf_h=NULL; /* RWF on CPU (2D) */
  rwf_h = sf_complexalloc( nx * nw );
  sf_complexread(rwf_h, nx*nw, Frwf);

  cuComplex *rwf_d; /* RWF on GPU (2D) */
  cudaMalloc((void **)&rwf_d, nx*nw*sizeof(cuComplex) );
  cudaMemcpy( rwf_d, rwf_h, nx*nw*sizeof(cuComplex), cudaMemcpyHostToDevice );
  
  /* Initialize Adjoint wavefield volumes on GPUs and set to zero */
  cuComplex *swfadj_d;
  cudaMalloc((void **)&swfadj_d, nx*nw*sizeof(cuComplex)); 
  cudaMemset(swfadj_d, 0, nx*nw*sizeof(cuComplex));

  cuComplex *rwfadj_d;
  cudaMalloc((void **)&rwfadj_d, nx*nw*sizeof(cuComplex)); 
  cudaMemset(rwfadj_d, 0, nx*nw*sizeof(cuComplex));

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
  float *vmin_h=NULL;
  vmin_h = sf_floatalloc(nz);
  float m;
  for (int  iz=0; iz<nz; iz++) {
	m = SF_HUGE;  
	for (int ix=0; ix<nx; ix++) {
		if (vel_h[iz*nx+ix] < m) m = vel_h[iz*nx+ix];
	}
	vmin_h[iz]=m;
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
  /* Handle tridiagonal solver. */
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
  for (int ixh=0; ixh < nh*nx; ixh++) {
  	xigslice_h[ixh] = xig_h[(nz-1)*nx*nh+ixh];
  }
  cudaMemcpy(xig_d, xigslice_h, nx*nh*sizeof(cuComplex), cudaMemcpyHostToDevice );

  /**************************************************************/
  // 
  // . . COMPUTE ADJOINT WAVEFIELDS
  //
  /**************************************************************/	
  make_adj_wflds<<<dimGrid,dimBlock>>>(xig_d,swf_d,rwf_d,swfadj_d,rwfadj_d,nh,nx,nw);

  /**************************************************************/
  // 
  // . . COMPUTE GRADIENT
  //
  /**************************************************************/	
  compute_gradient<<<dimGrid2,dimBlock2>>>(grd_d,swf_d,rwf_d,swfadj_d,rwfadj_d,nx,ow,dw,nw);		

  /**************************************************************/
  // 
  // . . COPY BACK GRADIENT TO HOST
  //
  /**************************************************************/	
  cudaMemcpy( grdslice_h, grd_d, nx*sizeof(float), cudaMemcpyDeviceToHost );
  for (int ix=0; ix<nx; ix++) {
  	grd_h[(nz-1)*nx+ix]=grdslice_h[ix];
  }
  		
  for (int iz = nz-2; iz > -1; iz--) {  /* Depth Loop */ 
    	if (verbose) fprintf(stderr,"\b\b\b\b\b\b\b\b\b\b\b\b%d of %d",nz-iz,nz); 

  		for (int ixh=0; ixh < nh*nx; ixh++) {
  			xigslice_h[ixh] = xig_h[iz*nx*nh+ixh];
  		}

   		/**************************************************************/
  		// 
  		// . . WAVEFIELD TAPER
  		//
  		/**************************************************************/
  		wfld_taper<<<dimGrid, dimBlock >>>(swf_d,   tap_d,nx,nw);
		wfld_taper<<<dimGrid, dimBlock >>>(rwf_d,   tap_d,nx,nw);
  		wfld_taper<<<dimGrid, dimBlock >>>(swfadj_d,tap_d,nx,nw);
		wfld_taper<<<dimGrid, dimBlock >>>(rwfadj_d,tap_d,nx,nw);

   		/**************************************************************/
  		// 
  		// . . SPLIT-STEP PROPAGATION
  		//
  		/**************************************************************/
  		prop_SSF<<<dimGrid,dimBlock>>>(swf_d,   vel_d,dw,ow,acaus,iz,dz,nx,nw);
		prop_SSF<<<dimGrid,dimBlock>>>(rwf_d,   vel_d,dw,ow, caus,iz,dz,nx,nw);
  		prop_SSF<<<dimGrid,dimBlock>>>(swfadj_d,vel_d,dw,ow,acaus,iz,dz,nx,nw);
		prop_SSF<<<dimGrid,dimBlock>>>(rwfadj_d,vel_d,dw,ow, caus,iz,dz,nx,nw);

  		/**************************************************************/
  		// 
  		// . . STEP WAVEFIELDS using tridiagonal and CUSPARSE
  		//
  		/**************************************************************/
		setup_FD<<<dimGrid,dimBlock>>>(swf_d,v_d,rax_d,rbx_d,rcx_d,vel_d,dw,ow,-acaus*aabb_h[0],aabb_h[1],nx,nw,iz);
#if (CUDART_VERSION >= 10000)
		cusparseCgtsv2StridedBatch_bufferSizeExt(cusparseHandle, nx, rcx_d, rax_d, rbx_d, v_d, nw, nx, &pbuffersize);
		cusparseCgtsv2StridedBatch(cusparseHandle, nx, rcx_d, rax_d, rbx_d, v_d, nw, nx, &pbuffersize);
#else
		cusparseCgtsvStridedBatch(cusparseHandle, nx, rcx_d, rax_d, rbx_d, v_d, nw, nx);
#endif
		copy_wfld<<<dimGrid,dimBlock>>>(v_d,swf_d,nx,nw); /* SWF */
		setup_FD<<<dimGrid,dimBlock>>>(swf_d,v_d,rax_d,rbx_d,rcx_d,vel_d,dw,ow,-acaus*aabb_h[2],aabb_h[3],nx,nw,iz);
#if (CUDART_VERSION >= 10000)
		cusparseCgtsv2StridedBatch_bufferSizeExt(cusparseHandle, nx, rcx_d, rax_d, rbx_d, v_d, nw, nx, &pbuffersize);
		cusparseCgtsv2StridedBatch(cusparseHandle, nx, rcx_d, rax_d, rbx_d, v_d, nw, nx, &pbuffersize);
#else
		cusparseCgtsvStridedBatch(cusparseHandle, nx, rcx_d, rax_d, rbx_d, v_d, nw, nx);
#endif
		copy_wfld<<<dimGrid,dimBlock>>>(v_d,swf_d,nx,nw);

		setup_FD<<<dimGrid,dimBlock>>>(rwf_d,v_d,rax_d,rbx_d,rcx_d,vel_d,dw,ow, -caus*aabb_h[0],aabb_h[1],nx,nw,iz);
#if (CUDART_VERSION >= 10000)
		cusparseCgtsv2StridedBatch_bufferSizeExt(cusparseHandle, nx, rcx_d, rax_d, rbx_d, v_d, nw, nx, &pbuffersize);
		cusparseCgtsv2StridedBatch(cusparseHandle, nx, rcx_d, rax_d, rbx_d, v_d, nw, nx, &pbuffersize);
#else
		cusparseCgtsvStridedBatch(cusparseHandle, nx, rcx_d, rax_d, rbx_d, v_d, nw, nx);
#endif
		copy_wfld<<<dimGrid,dimBlock>>>(v_d,rwf_d,nx,nw); /* RWF */
		setup_FD<<<dimGrid,dimBlock>>>(rwf_d,v_d,rax_d,rbx_d,rcx_d,vel_d,dw,ow, -caus*aabb_h[2],aabb_h[3],nx,nw,iz);
#if (CUDART_VERSION >= 10000)
		cusparseCgtsv2StridedBatch_bufferSizeExt(cusparseHandle, nx, rcx_d, rax_d, rbx_d, v_d, nw, nx, &pbuffersize);
		cusparseCgtsv2StridedBatch(cusparseHandle, nx, rcx_d, rax_d, rbx_d, v_d, nw, nx, &pbuffersize);
#else
		cusparseCgtsvStridedBatch(cusparseHandle, nx, rcx_d, rax_d, rbx_d, v_d, nw, nx);
#endif
		copy_wfld<<<dimGrid,dimBlock>>>(v_d,rwf_d,nx,nw);
 
		setup_FD<<<dimGrid,dimBlock>>>(swfadj_d,v_d,rax_d,rbx_d,rcx_d,vel_d,dw,ow,-acaus*aabb_h[0],aabb_h[1],nx,nw,iz);
#if (CUDART_VERSION >= 10000)
		cusparseCgtsv2StridedBatch_bufferSizeExt(cusparseHandle, nx, rcx_d, rax_d, rbx_d, v_d, nw, nx, &pbuffersize);
		cusparseCgtsv2StridedBatch(cusparseHandle, nx, rcx_d, rax_d, rbx_d, v_d, nw, nx, &pbuffersize);
#else
		cusparseCgtsvStridedBatch(cusparseHandle, nx, rcx_d, rax_d, rbx_d, v_d, nw, nx);
#endif
		copy_wfld<<<dimGrid,dimBlock>>>(v_d,swfadj_d,nx,nw); /* SWFADJ */
		setup_FD<<<dimGrid,dimBlock>>>(swfadj_d,v_d,rax_d,rbx_d,rcx_d,vel_d,dw,ow,-acaus*aabb_h[2],aabb_h[3],nx,nw,iz);
#if (CUDART_VERSION >= 10000)
		cusparseCgtsv2StridedBatch_bufferSizeExt(cusparseHandle, nx, rcx_d, rax_d, rbx_d, v_d, nw, nx, &pbuffersize);
		cusparseCgtsv2StridedBatch(cusparseHandle, nx, rcx_d, rax_d, rbx_d, v_d, nw, nx, &pbuffersize);
#else
		cusparseCgtsvStridedBatch(cusparseHandle, nx, rcx_d, rax_d, rbx_d, v_d, nw, nx);
#endif
		copy_wfld<<<dimGrid,dimBlock>>>(v_d,swfadj_d,nx,nw);

		setup_FD<<<dimGrid,dimBlock>>>(rwfadj_d,v_d,rax_d,rbx_d,rcx_d,vel_d,dw,ow, -caus*aabb_h[0],aabb_h[1],nx,nw,iz);
#if (CUDART_VERSION >= 10000)
		cusparseCgtsv2StridedBatch_bufferSizeExt(cusparseHandle, nx, rcx_d, rax_d, rbx_d, v_d, nw, nx, &pbuffersize);
		cusparseCgtsv2StridedBatch(cusparseHandle, nx, rcx_d, rax_d, rbx_d, v_d, nw, nx, &pbuffersize);
#else
		cusparseCgtsvStridedBatch(cusparseHandle, nx, rcx_d, rax_d, rbx_d, v_d, nw, nx);
#endif
		copy_wfld<<<dimGrid,dimBlock>>>(v_d,rwfadj_d,nx,nw); /* RWFADJ */
		setup_FD<<<dimGrid,dimBlock>>>(rwfadj_d,v_d,rax_d,rbx_d,rcx_d,vel_d,dw,ow, -caus*aabb_h[2],aabb_h[3],nx,nw,iz);
#if (CUDART_VERSION >= 10000)
		cusparseCgtsv2StridedBatch_bufferSizeExt(cusparseHandle, nx, rcx_d, rax_d, rbx_d, v_d, nw, nx, &pbuffersize);
		cusparseCgtsv2StridedBatch(cusparseHandle, nx, rcx_d, rax_d, rbx_d, v_d, nw, nx, &pbuffersize);
#else
		cusparseCgtsvStridedBatch(cusparseHandle, nx, rcx_d, rax_d, rbx_d, v_d, nw, nx);
#endif
		copy_wfld<<<dimGrid,dimBlock>>>(v_d,rwfadj_d,nx,nw);
  
   		/**************************************************************/
  		// 
  		// . . FORWARD FOURIER TRANSFORM
  		//
  		/**************************************************************/
 		cufftExecC2C(plan, swf_d,    swf_d,    CUFFT_FORWARD);
		cufftExecC2C(plan, rwf_d,    rwf_d,    CUFFT_FORWARD);
		cufftExecC2C(plan, swfadj_d, swfadj_d, CUFFT_FORWARD);
		cufftExecC2C(plan, rwfadj_d, rwfadj_d, CUFFT_FORWARD);
	
   		/**************************************************************/
  		// 
  		// . . PHASE CORRECTIONS
  		//
  		/**************************************************************/
 		phase_correction<<<dimGrid, dimBlock>>>(swf_d,   vmin_h[iz], dw, ow, kxmap_d, iz,-acaus, nx, nw, dz);
		phase_correction<<<dimGrid, dimBlock>>>(rwf_d,   vmin_h[iz], dw, ow, kxmap_d, iz,-  caus, nx, nw, dz);
 		phase_correction<<<dimGrid, dimBlock>>>(swfadj_d,vmin_h[iz], dw, ow, kxmap_d, iz,- acaus, nx, nw, dz);
		phase_correction<<<dimGrid, dimBlock>>>(rwfadj_d,vmin_h[iz], dw, ow, kxmap_d, iz,-  caus, nx, nw, dz);
	
   		/**************************************************************/
  		// 
  		// . . INVERSE FOURIER TRANSFORM
  		//
  		/**************************************************************/
  		cufftExecC2C(plan, swf_d,    swf_d,    CUFFT_INVERSE);
		cufftExecC2C(plan, rwf_d,    rwf_d,    CUFFT_INVERSE);
  		cufftExecC2C(plan, swfadj_d, swfadj_d, CUFFT_INVERSE);
		cufftExecC2C(plan, rwfadj_d, rwfadj_d, CUFFT_INVERSE);

   		/**************************************************************/
  		// 
  		// . . RESCALE WAVEFIELDS
  		//
  		/**************************************************************/
  		rescale_wfld<<< dimGrid, dimBlock >>>(swf_d,   nx,nw);
		rescale_wfld<<< dimGrid, dimBlock >>>(rwf_d,   nx,nw);
  		rescale_wfld<<< dimGrid, dimBlock >>>(swfadj_d,nx,nw);
		rescale_wfld<<< dimGrid, dimBlock >>>(rwfadj_d,nx,nw);
	
   		/**************************************************************/
  		// 
  		// . . WAVEFIELD TAPER
  		//
  		/**************************************************************/
  		wfld_taper<<<dimGrid, dimBlock>>>(swf_d,   tap_d,nx,nw);
		wfld_taper<<<dimGrid, dimBlock>>>(rwf_d,   tap_d,nx,nw);
  		wfld_taper<<<dimGrid, dimBlock>>>(swfadj_d,tap_d,nx,nw);
		wfld_taper<<<dimGrid, dimBlock>>>(rwfadj_d,tap_d,nx,nw);
	
   		/**************************************************************/
  		// 
  		// . . COMPUTE ADJOINT WAVEFIELDS
  		//
  		/**************************************************************/	

  		cudaMemcpy(xig_d, xigslice_h, nx*nh*sizeof(float), cudaMemcpyHostToDevice );
  		
  		/* Make the adjoint wavefields */
  		make_adj_wflds<<<dimGrid, dimBlock>>>(xig_d,swf_d,rwf_d,swfadj_d,rwfadj_d,nh,nx,nw);
	
   		/**************************************************************/
  		// 
  		// . . COMPUTE GRADIENT 
  		//
  		/**************************************************************/	
  		compute_gradient<<<dimGrid2,dimBlock2>>>(grd_d,swf_d,rwf_d,swfadj_d,rwfadj_d,nx,ow,dw,nw);		

		/* Copy back XIG slice and write output (include sum over frequency) */
  		cudaMemcpy( grdslice_h, grd_d, nx*sizeof(float), cudaMemcpyDeviceToHost );
  		for (int ix=0; ix<nx; ix++) {
  			grd_h[iz*nx+ix]=grdslice_h[ix];
  		}	
  		
    } /* END DEPTH LOOP */

  /* Write out gradient */
  if (verbose) fprintf(stderr,"\n"); 
  sf_warning("WRITING OUTPUT GRADIENT");
  
  sf_floatwrite(grd_h, nx*nz, Fgrd);

  /**************************************************************/
  // 
  // . . COMPUTE GRADIENT 
  //
  /**************************************************************/
  cufftDestroy(plan);
  cudaFree(vel_d);    cudaFree(xig_d); 
  cudaFree(grd_d);    cudaFree(tap_d);   cudaFree(kxmap_d);
  cudaFree(swf_d);    cudaFree(rwf_d);
  cudaFree(swfadj_d); cudaFree(rwfadj_d);
  cudaFree(rax_d); cudaFree(rbx_d); cudaFree(rcx_d);
  cudaFree(v_d);
  

  exit(0);
}
