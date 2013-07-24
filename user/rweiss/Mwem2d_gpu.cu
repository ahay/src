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
#include "cusparse_v2.h"
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
#define BLOCK1D 256
#define BLOCK2D 16

// checks the current GPU device for an error flag and prints to stderr
//static void sf_check_gpu_error (const char *msg) {
//    cudaError_t err = cudaGetLastError ();
//    if (cudaSuccess != err)
//        sf_error ("Cuda error: %s: %s", msg, cudaGetErrorString (err));
//}

// ENTRY POINT
int main(int argc, char* argv[])
{
  int nw,nx,nz,nh;
  float dx,dz,ow,dw,oh,dh;
  bool wantwf,wantillum;
  int nxtap;
  bool verbose;
  /* I/O files */
  sf_file Fvel  = NULL; /* velocity file */
  sf_file Fxig  = NULL; /* input XIG file */
  sf_file Fswf  = NULL; /* input SWF at iz=0 */
  sf_file Frwf  = NULL; /* input RWF at iz=0 */
  sf_file Fxigo = NULL; /* output xig file */
  sf_file Fswfo = NULL; /* output SWF at iz=NZ-1 file*/
  sf_file Frwfo = NULL; /* output RWF at iz=NZ-1 file*/
  sf_file Fsillum=NULL; /* output SWF illumination */
  sf_file Frillum=NULL; /* output RWF illumination */
  
  /*------------------------------------------------------------*/
  sf_init(argc,argv); /* init RSF */
  sf_axis ax,aw,az,ah,anull; /* cube axes */
  
  /* Set up file objects */
  Fvel = sf_input("vel"); /* input velocity */
  Fxig = sf_input("in");  sf_settype(Fxig,SF_COMPLEX);/* input xig */
  Fswf = sf_input("swf"); sf_settype(Fswf,SF_COMPLEX);/* INPUT SWF at iz=0 */
  Frwf = sf_input("rwf"); sf_settype(Frwf,SF_COMPLEX);/* INPUT RWF at iz=0 */
  Fxigo = sf_output("out");  sf_settype(Fxigo,SF_COMPLEX);/* OUTPUT XIG */
  Fswfo = sf_output("swfout"); sf_settype(Fswfo,SF_COMPLEX);/* OUTPUT SWF at iz=NZ-1 */
  Frwfo = sf_output("rwfout"); sf_settype(Frwfo,SF_COMPLEX);/* OUTPUT rWF at iz=NZ-1 */
  Fsillum=sf_output("sillum"); sf_settype(Fsillum,SF_COMPLEX);/* OUTPUT SWF ILLUMINATION */
  Frillum=sf_output("rillum"); sf_settype(Frillum,SF_COMPLEX);/* OUTPUT SWF ILLUMINATION */
  
  /*------------------------------------------------------------*/
  /* init GPU */
  //  int gpu;
  //if (! sf_getint("gpu", &gpu)) gpu = 0;	/* ID of the GPU to be used */
  //sf_warning("using GPU #%d", gpu);
  //cudaSetDevice(gpu);
  cudaSetDevice(0);
  //  cudaDeviceSetCacheConfig(cudaFuncCachePreferL1);
  sf_warning("FINISHED SETTING GPU");

   /* Read in parameter */ 
  if(! sf_getint("nxtap",&nxtap)) nxtap=40;  /* TAPER size */
  if(! sf_getbool("verbose",&verbose)) verbose=false; /* VERBOSITY flag */
  if(! sf_getbool("wantwf",&wantwf)) wantwf=false; /* Want output wavefields */
  if(! sf_getbool("wantillum",&wantillum)) wantillum=false; /* Want output wavefields */

  sf_warning("FINISHED PARAMS");

  /* Get axes information */
  ax = sf_iaxa(Frwf,1);
  aw = sf_iaxa(Frwf,2);
  az = sf_iaxa(Fvel,2);
  ah = sf_iaxa(Frwf,1);
  anull=sf_iaxa(Frwf,1);

  /* Pull out individual values */
  nx = sf_n(ax); dx = sf_d(ax);
  nw = sf_n(aw); dw = sf_d(aw); ow = sf_o(aw); 
  nz = sf_n(az); dz = sf_d(az);

  /* xig shift parameter */
  if(! sf_getint("nh",&nh)) nh=0;
  dh = dx;
  oh = -(float)nh*dh;
  if (nh==0) oh = 0.f;
  nh = 2*nh+1;

  /* Set up xig shift axis */
  sf_setn(ah,nh);
  sf_seto(ah,oh);
  sf_setd(ah,dh);
  sf_setn(anull,1);
  sf_setd(anull,1.);
  sf_seto(anull,0.);

  /* Output axes information */
  if (wantwf) {
    sf_oaxa(Fswfo,ax,1);
    sf_oaxa(Fswfo,aw,2);
    sf_oaxa(Fswfo,anull,3);

    sf_oaxa(Frwfo,ax,1);
    sf_oaxa(Frwfo,aw,2);
    sf_oaxa(Frwfo,anull,3);
  }

  if (wantillum) {
  	sf_oaxa(Fsillum,ax,1);
  	sf_oaxa(Fsillum,az,2);
  	sf_oaxa(Fsillum,anull,3);

  	sf_oaxa(Frillum,ax,1);
  	sf_oaxa(Frillum,az,2);
   	sf_oaxa(Frillum,anull,3);
 }
  
  /* Output axes information*/
  if (verbose) {
    sf_raxa(az); sf_raxa(ax);
    sf_raxa(aw); sf_raxa(ah);
  }

  sf_warning("FINISHED LOGIC");

  /**************************************************************/
  /* Read in velocity to CPU and copy to GPU */
  float *slo_h=NULL; /* s array (CPU) */
  float *vel_h=NULL; /* v array (CPU) */
  float *vel_d; /* v array (GPU) */ 
  slo_h = sf_floatalloc(nx*nz);
  vel_h = sf_floatalloc(nx*nz);
  sf_floatread( slo_h, nx*nz, Fvel );	
  
  // Convert from slowness to velocity
  for (int ixz=0; ixz < nz*nx; ixz++) {
    vel_h[ixz] = 1.f/slo_h[ixz];
  }

  cudaMalloc((void **)&vel_d, nx*nz*sizeof(float) );
  cudaMemcpy(vel_d, vel_h, nx*nz*sizeof(float), cudaMemcpyHostToDevice );
  sf_warning("FINISHED VELOCITY");
    
  /**************************************************************/
  /* Read in xig to CPU and copy to GPU */
  sf_complex *xig_h=NULL; /* XIG (CPU) 2D */
  cuComplex *xig_d; /* XIG array (GPU) - 2D slice*/
  cudaMalloc((void **)&xig_d, nx*nh*sizeof(cuComplex) );
  xig_h = sf_complexalloc( nx * nh );
  cudaMemset(xig_d,   0., nx*nh*sizeof(cuComplex));
  sf_warning("FINISHED XIG");

  /**************************************************************/
  /* Read in SWF to CPU and copy to GPU */	
  sf_complex *swf_h=NULL; /* SWF on CPU (2D) */
  cuComplex *swf_d; /* SWF on GPU (2D) */
  cudaMalloc((void **)&swf_d, nx*nw*sizeof(cuComplex) );
  swf_h = sf_complexalloc( nx * nw );
  sf_complexread(swf_h, nx*nw, Fswf);
  cudaMemcpy( swf_d, swf_h, nx*nw*sizeof(cuComplex), cudaMemcpyHostToDevice );
  sf_warning("FINISHED SWF");

  /**************************************************************/
  /* Read in RWF to CPU and copy to GPU */	
  sf_complex *rwf_h=NULL; /* RWF on CPU (2D) */
  cuComplex *rwf_d; /* RWF on GPU (2D) */
  cudaMalloc((void **)&rwf_d, nx*nw*sizeof(cuComplex) );
  rwf_h = sf_complexalloc( nx * nw );
  sf_complexread(rwf_h, nx*nw, Frwf);
  cudaMemcpy( rwf_d, rwf_h, nx*nw*sizeof(cuComplex), cudaMemcpyHostToDevice );
  sf_warning("FINISHED RWF");

  /**************************************************************/
  /* Initialize wavefield taper and copy to GPU */
  sf_complex *tap_h=NULL; /* taper on CPU */
  cuComplex *tap_d; /* taper on GPU */
  tap_h = sf_complexalloc(nx);
  int j1;
  for (int ix=nxtap; ix<nx; ix++) {
    tap_h[ix]= sf_cmplx(1.f,0.f);
  }
  for (int ix=0; ix < nxtap; ix++) {
    j1 = abs(nxtap-ix-1.f);
    tap_h[ix] = sf_cmplx( cos ( SF_PI/2.f*(float)j1/((float)nxtap-1.f) ), 0.f);
  }
  
  for (int ix=0; ix < nxtap; ix++) {
    j1 = abs(nxtap-ix-1.f);
    tap_h[nx-ix-1] = sf_cmplx( cos (SF_PI/2.f*(float)j1/((float)nxtap-1.f) ),0.f);
  }
    
  tap_h[0]=sf_cmplx(0.f,0.f); 
  tap_h[nx-1]=sf_cmplx(0.f,0.f);
  
  cudaMalloc((void **)&tap_d, nx*sizeof(cuComplex) );
  cudaMemcpy( tap_d, tap_h, nx*sizeof(cuComplex), cudaMemcpyHostToDevice );
  sf_warning("FINISHED TAPER");

  /**************************************************************/
  /* Set up vmin array */
  float *vmin_h=NULL; 
  vmin_h = sf_floatalloc(nz);
  float m;
  int iz,ix;
  for (iz=0; iz<nz; iz++) {
    m = SF_HUGE;  
    for (ix=0; ix<nx; ix++) {
      m = SF_MIN(m,vel_h[iz*nx+ix]);
    }
    vmin_h[iz]=m;
  }	
  sf_warning("FINISHED MIN VEL");
  
  /**************************************************************/
  /* Set up illumination arrays */
  sf_complex *Sillum_h=NULL;  
  sf_complex *Rillum_h=NULL; 
  cuComplex *illum_d;
  cudaMalloc((void **)&illum_d, nx*sizeof(cuComplex) );
  Sillum_h = sf_complexalloc(nx);
  Rillum_h = sf_complexalloc(nx);
  sf_warning("FINISHED MIN VEL");  
  
  /**************************************************************/
  /* Set up Fourier domain wavenumbers */	
  float *kxmap_h=NULL; float *kxmap_d;
  kxmap_h = sf_floatalloc(nx); 
  int mm=nx/2;
  float kscale = 2.f*SF_PI/( (float)nx * dx);
  int ii;
  for (ii=0; ii<mm+1; ii++) {
    kxmap_h[ii]=kscale*((float)ii);
  }	

  for (ii=mm; ii<nx; ii++) {
    kxmap_h[ii]=kscale*((float)ii-(float)nx);
  }
  cudaMalloc((void **)&kxmap_d, nx*sizeof(float) );
  cudaMemcpy( kxmap_d, kxmap_h, nx*sizeof(float), cudaMemcpyHostToDevice );
  sf_warning("FINISHED WAVENUMBER");

 
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
  dim3 dimGrid( ceil(nx/16.f), ceil(nw/16.f), 1 ); 
  dim3 dimBlock2( 256, 1, 1);
  dim3 dimGrid2( ceil(nx/256.f), 1, 1);

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
  cudaMemset(rax_d, 0., nx*nw*sizeof(cuComplex));
  cudaMemset(rbx_d, 0., nx*nw*sizeof(cuComplex));
  cudaMemset(rcx_d, 0., nx*nw*sizeof(cuComplex));
  cudaMemset(v_d,   0., nx*nw*sizeof(cuComplex));
  sf_warning("FINISHED TEMP ALLOCATION");

  /**************************************************************/
  /* FFT set up Create a 1D FFT plan. */
  cufftHandle plan;
  cufftPlan1d(&plan, nx, CUFFT_C2C, nw);
  sf_warning("FINISHED FFT PLAN");

  /**************************************************************/
  /* Handle tridiagonal solver. */
  cusparseHandle_t cusparseHandle = 0;
  /* cusparseStatus_t cusparseStatus;
     cusparseStatus = cusparseCreate(&cusparseHandle); */
  sf_warning("FINISHED TRIDIAG SOLVER PLAN");

  /**************************************************************/
  /* MAIN LOOP */	
  for (int iz = 0; iz < nz; iz++) {  /* Depth Loop */ 
    if (verbose) fprintf(stderr,"\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%d of %d  ",iz+1,nz); 
        
    /* Split-step bulk-phase shift */
    prop_SSF<<<dimGrid,dimBlock>>>(swf_d,vel_d,dw,ow, caus,iz,dz,nx,nw);
    prop_SSF<<<dimGrid,dimBlock>>>(rwf_d,vel_d,dw,ow,acaus,iz,dz,nx,nw);
    
    /* Set up FD coefficients for SWF and solve using CUSPARSE*/
    /* FIRST STEP*/
    setup_FD<<<dimGrid,dimBlock>>>(swf_d,v_d,rax_d,rbx_d,rcx_d,vel_d,dw,ow,-caus*aabb_h[0],aabb_h[1],nx,nw,iz);
    cusparseCgtsvStridedBatch(cusparseHandle, nx, rcx_d, rax_d, rbx_d, v_d, nw, nx);
    copy_wfld<<<dimGrid,dimBlock>>>(v_d,swf_d,nx,nw);
    
    /* SECOND STEP*/
    setup_FD<<<dimGrid,dimBlock>>>(swf_d,v_d,rax_d,rbx_d,rcx_d,vel_d,dw,ow,-caus*aabb_h[2],aabb_h[3],nx,nw,iz);
    cusparseCgtsvStridedBatch(cusparseHandle, nx, rcx_d, rax_d, rbx_d, v_d, nw, nx);
    copy_wfld<<<dimGrid,dimBlock>>>(v_d,swf_d,nx,nw);
    
    /* Set up FD coefficients for RWF and solve using CUSPARSE*/
    /* FIRST STEP*/
    setup_FD<<<dimGrid,dimBlock>>>(rwf_d,v_d,rax_d,rbx_d,rcx_d,vel_d,dw,ow,-acaus*aabb_h[0],aabb_h[1],nx,nw,iz);
    cusparseCgtsvStridedBatch(cusparseHandle, nx, rcx_d, rax_d, rbx_d, v_d, nw, nx);
    copy_wfld<<<dimGrid,dimBlock>>>(v_d,rwf_d,nx,nw);
    
    /* SECOND STEP*/
    setup_FD<<<dimGrid,dimBlock>>>(rwf_d,v_d,rax_d,rbx_d,rcx_d,vel_d,dw,ow,-acaus*aabb_h[2],aabb_h[3],nx,nw,iz);
    cusparseCgtsvStridedBatch(cusparseHandle, nx, rcx_d, rax_d, rbx_d, v_d, nw, nx);
    copy_wfld<<<dimGrid,dimBlock>>>(v_d,rwf_d,nx,nw);
    
    /* High-angle filter FFT to kx */
    cufftExecC2C(plan, swf_d, swf_d, CUFFT_FORWARD);
    cufftExecC2C(plan, rwf_d, rwf_d, CUFFT_FORWARD);
    
    /* Fourier domain filtering */
    phase_correction<<<dimGrid, dimBlock>>>(swf_d, vmin_h[iz], dw, ow, kxmap_d, iz, acaus, nx, nw, dz);
    phase_correction<<<dimGrid, dimBlock>>>(rwf_d, vmin_h[iz], dw, ow, kxmap_d, iz,  caus, nx, nw, dz);
    
    /* High-angle filter FFT to kx */
    cufftExecC2C(plan, swf_d, swf_d, CUFFT_INVERSE);
    cufftExecC2C(plan, rwf_d, rwf_d, CUFFT_INVERSE);
    
    /* Rescale wavefields by nx */
    rescale_wfld<<< dimGrid, dimBlock >>>(swf_d,nx,nw);
    rescale_wfld<<< dimGrid, dimBlock >>>(rwf_d,nx,nw);
    
    /* WAVEFIELD TAPER */
    wfld_taper<<<dimGrid, dimBlock>>>(swf_d,tap_d,nx,nw);
    wfld_taper<<<dimGrid, dimBlock>>>(rwf_d,tap_d,nx,nw);
    
    /* XIG Imaging Condition */
    if (nh==1) {  
      imaging_condition<<<dimGrid2,dimBlock2>>>(swf_d,rwf_d,xig_d,iz,nh,nx,nw);
    } else {
      xig_condition<<<dimGrid2,dimBlock2>>>(swf_d,rwf_d,xig_d,iz,nh,nx,nw);
    }
    
	/* Copy back XIG slice and write output (include sum over frequency) */
	cudaMemcpy( xig_h, xig_d, nx*nh*sizeof(cuComplex), cudaMemcpyDeviceToHost );
	sf_complexwrite(xig_h, nx*nh, Fxigo);
    
    if (wantillum) {
	    /* Compute Source Illumination */
	    illumination<<<dimGrid2,dimBlock2>>>(swf_d,illum_d,iz,nx,nw);
	    cudaMemcpy( Sillum_h, illum_d, nx*sizeof(cuComplex), cudaMemcpyDeviceToHost );
    	sf_complexwrite( Sillum_h, nx, Fsillum );

	    /* Compute Receiver Illumination */
	    illumination<<<dimGrid2,dimBlock2>>>(rwf_d,illum_d,iz,nx,nw);
	    cudaMemcpy( Rillum_h, illum_d, nx*sizeof(cuComplex), cudaMemcpyDeviceToHost );
    	sf_complexwrite( Rillum_h, nx, Frillum );
	}
    
  } /* END DEPTH */
  
  if (wantwf) {
    /* Copy back SWF and RWF and write output */
    sf_warning("WRITING OUT RWF and SWF");
    cudaMemcpy( swf_h, swf_d, nx*nw*sizeof(cuComplex), cudaMemcpyDeviceToHost );
    sf_complexwrite( swf_h, nx*nw, Fswfo);
    cudaMemcpy( rwf_h, rwf_d, nx*nw*sizeof(cuComplex), cudaMemcpyDeviceToHost );
    sf_complexwrite( rwf_h, nx*nw, Frwfo);
  }
  
  cufftDestroy(plan);
  cudaFree(vel_d); cudaFree(xig_d); cudaFree(swf_d); cudaFree(rwf_d);
  cudaFree(tap_d); cudaFree(kxmap_d); cudaFree(illum_d);
  
  exit(0);
}
