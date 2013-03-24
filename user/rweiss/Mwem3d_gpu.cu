/* 3D ISOTROPIC wave-equation finite-difference migration on GPU */
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
//#include <cuda_runtime_api.h>
#include "cusparse_v2.h"
#include <cuComplex.h>
#include <cufft.h>

extern "C" {
#include <rsf.h>
}

#include "wem3d_kernels.cu"
#define a1 0.040315157
#define b1 0.873981642
#define a2 0.457289566
#define b2 0.222691983
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define BLOCK2D 16
#define DNHX 1
#define DNHY 1

// checks the current GPU device for an error flag and prints to stderr
/*
static void sf_check_gpu_error (const char *msg) {
    cudaError_t err = cudaGetLastError ();
    if (cudaSuccess != err)
        sf_error ("Cuda error: %s: %s", msg, cudaGetErrorString (err));
}
*/

// ENTRY POINT
int main(int argc, char* argv[])
{
  int nw,nx,nz,ny;
  float dx,dz,dy,ow,dw,ohy,ohx,dhx,dhy;
  bool wantwf;
  int nxtap,nytap;
  bool verbose;
  /* I/O files */
  sf_file Fvel  = NULL; /* velocity file */
  sf_file Fxig  = NULL; /* input XIG file */
  sf_file Fswf  = NULL; /* input SWF at iz=0 */
  sf_file Frwf  = NULL; /* input RWF at iz=0 */
  sf_file Fxigo = NULL; /* output xig file */
  sf_file Fswfo = NULL; /* output SWF at iz=NZ-1 file*/
  sf_file Frwfo = NULL; /* output RWF at iz=NZ-1 file*/
  sf_file Fkxmap = NULL; /* output RWF at iz=NZ-1 file*/
 
  /*------------------------------------------------------------*/
  sf_init(argc,argv); /* init RSF */
  sf_axis ax,ay,aw,az,ahx,ahy,anull; /* cube axes */
  
  /* Set up file objects */
  Fxig = sf_input("in"); /* input xig */
  Fvel = sf_input("vel"); /* input velocity */
  Fswf = sf_input("swf"); sf_settype(Fswf,SF_COMPLEX);/* INPUT SWF at iz=0 */
  Frwf = sf_input("rwf"); sf_settype(Frwf,SF_COMPLEX);/* INPUT RWF at iz=0 */
  Fswfo = sf_output("swfout"); sf_settype(Fswfo,SF_COMPLEX);/* OUTPUT SWF at iz=NZ-1 */
  Frwfo = sf_output("rwfout"); sf_settype(Frwfo,SF_COMPLEX);/* OUTPUT rWF at iz=NZ-1 */
  Fkxmap = sf_output("kxmap"); /* OUTPUT rWF at iz=NZ-1 */
  Fxigo = sf_output("out"); /* OUTPUT XIG */

  /*------------------------------------------------------------*/
  /* init GPU */
  int gpu;
  if (! sf_getint("gpu", &gpu)) gpu = 0;	/* ID of the GPU to be used */
  sf_warning("using GPU #%d", gpu);
  cudaSetDevice(gpu);
  cudaDeviceSetCacheConfig(cudaFuncCachePreferL1);

   /* Read in parameter */ 
  if(! sf_getint("nxtap",&nxtap)) nxtap=40;  /* TAPER size */
  if(! sf_getint("nytap",&nytap)) nytap=40;  /* TAPER size */
  if(! sf_getbool("verbose",&verbose)) verbose=false; /* VERBOSITY flag */
  if(! sf_getbool("wantwf",&wantwf)) wantwf=false; /* Want output wavefields */

  /* Get axes information */
  ax = sf_iaxa(Frwf,1);
  ay = sf_iaxa(Frwf,2);
  aw = sf_iaxa(Frwf,3);
  az = sf_iaxa(Fvel,3);
  ahx = sf_iaxa(Frwf,1);
  ahy = sf_iaxa(Frwf,2);
  anull=sf_iaxa(Frwf,1);

  /* Pull out individual values */
  nx = sf_n(ax); dx = sf_d(ax);
  ny = sf_n(ay); dy = sf_d(ay);
  nw = sf_n(aw); dw = sf_d(aw); ow = sf_o(aw); 
  nz = sf_n(az); dz = sf_d(az);

  /* xig shift parameter */
  dhx = dx;  ohx = -(float)(DNHX-1)*dhx; 
  dhy = dy;  ohy = -(float)(DNHY-1)*dhy; 

  /* Set up xig shift axis */
  sf_setn(ahx,DNHY);  sf_seto(ahx,ohx);  sf_setd(ahx,dhx);
  sf_setn(ahy,DNHY);  sf_seto(ahy,ohy);  sf_setd(ahy,dhy);    
  sf_setn(anull,1);  sf_setd(anull,1.);  sf_seto(anull,0.);

  sf_oaxa(Fkxmap,ax,1);
  sf_oaxa(Fkxmap,anull,2);
  sf_oaxa(Fkxmap,anull,3);
  sf_oaxa(Fkxmap,anull,4);
  sf_oaxa(Fkxmap,anull,5);

  /* Output axes information */
  if (wantwf) {
    sf_oaxa(Fswfo,ax,1);
    sf_oaxa(Fswfo,ay,2);
    sf_oaxa(Fswfo,aw,3);
    sf_oaxa(Fswfo,anull,4);
    sf_oaxa(Fswfo,anull,5);
    sf_oaxa(Frwfo,ax,1);
    sf_oaxa(Frwfo,ay,2);
    sf_oaxa(Frwfo,aw,3);
    sf_oaxa(Frwfo,anull,4);
    sf_oaxa(Frwfo,anull,5);
  }

  /* Output axes information*/
  if (verbose) {
    sf_raxa(aw); 
    sf_raxa(ahx);
    sf_raxa(ahy);
  	sf_raxa(ax);
  	sf_raxa(ay);
  	sf_raxa(az); 
  }

  /**************************************************************/
  /* Read in velocity to CPU and copy to GPU */
  float *vel_h=NULL; /* v array (CPU) */
  float *vel_d; /* v array (GPU) */ 
  vel_h = sf_floatalloc(nx*ny);
  cudaMalloc((void **)&vel_d, ny*nx*sizeof(float) );

  /**************************************************************/
  /* Read in xig to CPU and copy to GPU */
  float *xig_h=NULL; /* XIG (CPU) 2D */
  float *xig_d; /* XIG array (GPU) - 2D slice*/
  xig_h = sf_floatalloc( nx * ny * DNHX * DNHY );
  cudaMalloc((void **)&xig_d, nx*ny*DNHX*DNHY*sizeof(float) );
  cudaMemcpy( xig_d, xig_h, nx*ny*DNHX*DNHY*sizeof(float), cudaMemcpyHostToDevice );

  /**************************************************************/
  /* Read in SWF to CPU and copy to GPU */	
  sf_complex *swf_h=NULL; /* SWF on CPU (3D) */
  cuComplex *swf_d; /* SWF on GPU (3DD) */
  swf_h = sf_complexalloc( nx*ny*nw );
  sf_complexread(swf_h, nx*ny*nw, Fswf);
  cudaMalloc((void **)&swf_d, nx*ny*nw*sizeof(cuComplex) );
  cudaMemcpy( swf_d, swf_h, nx*ny*nw*sizeof(cuComplex), cudaMemcpyHostToDevice );

  /**************************************************************/
  /* Read in SWF to CPU and copy to GPU */	
  sf_complex *rwf_h=NULL; /* RWF on CPU (2D) */
  cuComplex *rwf_d; /* RWF on GPU (2D) */
  rwf_h = sf_complexalloc( nx * ny * nw );  
  sf_complexread(rwf_h, nx*ny*nw, Frwf);
  cudaMalloc((void **)&rwf_d, nx*ny*nw*sizeof(cuComplex) );
  cudaMemcpy( rwf_d, rwf_h, nx*ny*nw*sizeof(cuComplex), cudaMemcpyHostToDevice );

  /**************************************************************/
  /* Initialize wavefield taper and copy to GPU */
  sf_complex **tap_h=NULL; /* taper on CPU */
  sf_complex *tap1d_h=NULL; /* taper on CPU */
  cuComplex *tap_d; /* taper on GPU */
  tap_h = sf_complexalloc2(nx,ny);
  tap1d_h = sf_complexalloc(nx*ny);
  int j1,j2;
  
  /* Set all to unity */
  for (int iy=0; iy < ny; iy++) {
  	for (int ix=0; ix < nx; ix++) {
  		tap_h[iy][ix]= sf_cmplx(1.f,0.f);
  	}
  }
  
  /* Do upper left corner */
  for (int iy=0; iy < nytap; iy++) {
  	j2 = abs(nytap-iy-1.f);
  	
  	float ycomp = cos ( SF_PI/2.f*(float)j2 / ((float)nytap-1.f) );
  	for (int ix=0; ix < nxtap; ix++) {
	  	j1 = abs(nxtap-ix-1.f);
	  	float xcomp = cos ( SF_PI/2.f*(float)j1/((float)nxtap-1.f) );
	    tap_h[iy][ix] = sf_cmplx( xcomp * ycomp, 0.f);
  	}
  }

  /* Do upper right corner */
  for (int iy=0; iy < nytap; iy++) {
  	j2 = abs(nytap-iy-1.f);
  	float ycomp = cos ( SF_PI/2.f*(float)j2/((float)nytap-1.f) );
  	for (int ix=0; ix < nxtap; ix++) {
	  	j1 = abs(nxtap-ix-1.f);
	  	float xcomp = cos ( SF_PI/2.f*(float)j1/((float)nxtap-1.f) );
	    tap_h[iy][nx-ix-1] = sf_cmplx( xcomp * ycomp , 0.f);
  	}
  }  
  
  /* Do bottom left corner */
  for (int iy=0; iy < nytap; iy++) {
  	j2 = abs(nytap-iy-1.f);

  	float ycomp = cos ( SF_PI/2.f*(float)j2 / ((float)nytap-1.f) );

  	for (int ix=0; ix < nxtap; ix++) {
	  	j1 = abs(nxtap-ix-1.f);
	  	float xcomp = cos ( SF_PI/2.f*(float)j1 / ((float)nxtap-1.f) );
	    tap_h[ny-iy-1][ix] = sf_cmplx( xcomp * ycomp , 0.f);
  	}
  }  

  /* Do bottom right corner */
  for (int iy=0; iy < nytap; iy++) {
  	j2 = abs(nytap-iy-1.f);
  	float ycomp = cos ( SF_PI/2.f*(float)j2/((float)nytap-1.f) );
  	for (int ix=0; ix < nxtap; ix++) {
	  	j1 = abs(nxtap-ix-1.f);
	  	float xcomp = cos ( SF_PI/2.f*(float)j1/((float)nxtap-1.f) );
	    tap_h[ny-iy-1][nx-ix-1] = sf_cmplx( xcomp * ycomp , 0.f);
  	}
  }  
    
  /* Do upper boundary */
  for (int iy=0; iy < nytap; iy++) {
  	j2 = abs(nytap-iy-1.f);
  	float ycomp = cos ( SF_PI/2.f*(float)j2/((float)nytap-1.f) );
  	for (int ix=nxtap; ix < nx-nxtap; ix++) {
	    tap_h[iy][ix] = sf_cmplx( ycomp , 0.f);
  	}
  }  
  
  /* Do lower boundary */
  for (int iy=0; iy < nytap; iy++) {
  	j2 = abs(nytap-iy-1.f);
  	float ycomp = cos ( SF_PI/2.f*(float)j2/((float)nytap-1.f) );
  	for (int ix=nxtap; ix < nx-nxtap; ix++) {
	    tap_h[ny-iy-1][ix] = sf_cmplx( ycomp , 0.f);
  	}
  }  

  /* Do left boundary */
  for (int iy=nytap; iy < ny-nytap; iy++) {
  	for (int ix=0; ix < nxtap; ix++) {
  		j1 = abs(nxtap-ix-1.f);
  		float xcomp = cos ( SF_PI/2.f*(float)j1/((float)nxtap-1.f) );
  		tap_h[iy][ix] = sf_cmplx( xcomp , 0.f);
  	}
  }  
    
  /* Do right boundary */
  for (int iy=nytap; iy < ny-nytap; iy++) {
  	for (int ix=0; ix < nxtap; ix++) {
  		j1 = abs(nxtap-ix-1.f);
  		float xcomp = cos ( SF_PI/2.f*(float)j1/((float)nxtap-1.f) );
  		tap_h[iy][nx-ix-1] = sf_cmplx( xcomp , 0.f);
  	}
  }   
  
  /* ensure zero along boundaries*/
  for (int ix=0; ix < nx; ix++) {
  	tap_h[0][ix] = sf_cmplx(0.f,0.f);    
  }
  for (int ix=0; ix < nx; ix++) {
  	tap_h[ny-1][ix] = sf_cmplx(0.f,0.f);    
  }
  for (int iy=0; iy < ny; iy++) {
  	tap_h[iy][0] = sf_cmplx(0.f,0.f);    
  }
  for (int iy=0; iy < ny; iy++) {
  	tap_h[iy][nx-1] = sf_cmplx(0.f,0.f);    
  }
  
  /* Copy to device */
  cudaMalloc((void **)&tap_d, nx*ny*sizeof(cuComplex) );

  /* Set all to unity */
  for (int iy=0; iy < ny; iy++) {
  	for (int ix=0; ix < nx; ix++) {
  		int addr = iy*nx+ix;
  		tap1d_h[addr]= tap_h[iy][ix];
  	}
  }
  cudaMemcpy( tap_d, tap1d_h, nx*ny*sizeof(cuComplex), cudaMemcpyHostToDevice );

  /**************************************************************/
  /* Set up Fourier domain wavenumbers */	
  float *kxmap_h=NULL; float *kxmap_d;
  kxmap_h = sf_floatalloc(nx); 
  int mmx=nx/2;
  float kxscale = 2.f*SF_PI/( (float)nx * dx);
  for (int ii=0; ii < mmx+1; ii++) {
  	kxmap_h[ii]=kxscale*((float)ii);
  }	
  for (int ii=mmx; ii < nx; ii++) {
  	kxmap_h[ii]=kxscale*((float)ii-(float)nx);
  }
  cudaMalloc((void **)&kxmap_d, nx*sizeof(float) );
  cudaMemcpy( kxmap_d, kxmap_h, nx*sizeof(float), cudaMemcpyHostToDevice );
 

  /**************************************************************/
  /* Set up Fourier domain wavenumbers */	
  float *kymap_h=NULL; float *kymap_d;
  kymap_h = sf_floatalloc(ny); 
  int mmy=ny/2;
  float kyscale = 2.f*SF_PI/( (float)ny * dy);
  for (int ii=0; ii < mmy+1; ii++) {
  	kymap_h[ii]=kyscale*((float)ii);
  }	
  for (int ii=mmy; ii < ny; ii++) {
  	kymap_h[ii]=kyscale*((float)ii-(float)ny);
  }
  cudaMalloc((void **)&kymap_d, ny*sizeof(float) );
  cudaMemcpy( kymap_d, kymap_h, ny*sizeof(float), cudaMemcpyHostToDevice );
 
  /**************************************************************/
  /* Set aabb array FD coeffs including dz and dx */
  float *ab_h=NULL;
  ab_h = sf_floatalloc(8);
  
  /* Coeffs in X direction */
  ab_h[0] = a1 * dz / (2.f*dx*dx);
  ab_h[1] = b1 / (dx * dx);
  ab_h[2] = a1 * dz / (2.f*dy*dy);
  ab_h[3] = b1 / (dy * dy);
  
  /* Coeffs in Y-direction */
  ab_h[4] = a2 * dz / (2.f*dx*dx);
  ab_h[5] = b2 / (dx * dx);
  ab_h[6] = a2 * dz / (2.f*dy*dy);
  ab_h[7] = b2 / (dy * dy);
  
  /**************************************************************/
  /* Set up Array calls */	
  dim3 dimBlock(16,16,1); 
  dim3 dimGrid(ceil(nx/16.f),ceil(ny/16.f),1); 
  
  dim3 dimBlock2(32,8,1); 
  dim3 dimGrid2(ceil(nx/32.f),ceil(ny/8.f),1); 
 
  
  /* Define causality */
  float  caus= 1.f;
  float acaus=-1.f;

  /**************************************************************/
  /* Set up arrays required by FD code */	
  cuComplex *SWFslice_d, *RWFslice_d;
  cuComplex *ra_d,*rb_d,*rc_d,*v_d;
  
  /* Slice to do */
  cudaMalloc((void **)&SWFslice_d,nx*ny*sizeof(cuComplex)); 
  cudaMemset(SWFslice_d,      0., nx*ny*sizeof(cuComplex));
  cudaMalloc((void **)&RWFslice_d,nx*ny*sizeof(cuComplex)); 
  cudaMemset(RWFslice_d,      0., nx*ny*sizeof(cuComplex));

  /* Allocate Memory on GPUs */
  cudaMalloc((void **)&ra_d,nx*ny*sizeof(cuComplex)); 
  cudaMalloc((void **)&rb_d,nx*ny*sizeof(cuComplex)); 
  cudaMalloc((void **)&rc_d,nx*ny*sizeof(cuComplex)); 
  cudaMalloc((void **)&v_d, nx*ny*sizeof(cuComplex)); 

  /* Set equal to zero */
  cudaMemset(ra_d, 0.,nx*ny*sizeof(cuComplex));
  cudaMemset(rb_d, 0.,nx*ny*sizeof(cuComplex));
  cudaMemset(rc_d, 0.,nx*ny*sizeof(cuComplex));
  cudaMemset(v_d , 0.,nx*ny*sizeof(cuComplex));

  /**************************************************************/
  /* FFT set up Create a 2D FFT plan. */
  cufftHandle plan;
  cufftPlan2d(&plan, nx, ny, CUFFT_C2C);

  /**************************************************************/
  /* Handle tridiagonal solver. */
  cusparseHandle_t cusparseHandleX = 0;
  cusparseStatus_t cusparseStatusX;
  cusparseStatusX = cusparseCreate(&cusparseHandleX);

  cusparseHandle_t cusparseHandleY = 0;
  cusparseStatus_t cusparseStatusY;
  cusparseStatusY = cusparseCreate(&cusparseHandleY);

  /**************************************************************/
  /* MAIN LOOP */	
  for (int iz = 0; iz < nz; iz++) {  /* Depth Loop */ 
    	if (verbose) fprintf(stderr,"\b\b\b\b\b\b\b\b\b\b\b\b%d of %d",iz+1,nz); 

		/* Copy over velocity at depth slice */
	  	sf_floatread( vel_h, nx*ny, Fvel );	
		cudaMemcpy(vel_d, vel_h, ny*nx*sizeof(float), cudaMemcpyHostToDevice );

	   	float minvel = SF_HUGE;  
		for (int ixy=0; ixy < nx*ny; ixy++) {
			if (vel_h[ixy] < minvel) minvel = vel_h[ixy];
		}

  		/* Zero XIG array */
	    cudaMemset(xig_d, 0.f, nx*ny*DNHX*DNHY*sizeof(float) );

  		for (int iw = 0; iw < nw; iw++) { /* Frequency Loop */
  		
			/* CURRENT FREQUENCY */  		
			float ww = 2.f*SF_PI*(((float)(iw-1))*dw+ow);

			/* Extract current wavefield slice */
			extract_slice<<<dimGrid, dimBlock>>>(swf_d,SWFslice_d,iw,nx,ny);
			extract_slice<<<dimGrid, dimBlock>>>(rwf_d,RWFslice_d,iw,nx,ny);

			/* WAVEFIELD TAPER */	
			wfld_taper3d<<<dimGrid, dimBlock >>>(SWFslice_d,tap_d,nx,ny);
			wfld_taper3d<<<dimGrid, dimBlock >>>(RWFslice_d,tap_d,nx,ny);

			/* Split-step bulk-phase shift */
			prop_SSF3d<<<dimGrid,dimBlock>>>(SWFslice_d,vel_d,ww, caus,dz,nx,ny);
			prop_SSF3d<<<dimGrid,dimBlock>>>(RWFslice_d,vel_d,ww,acaus,dz,nx,ny);

			/* Set up FD coefficients for SWF and solve using CUSPARSE*/
			for (int is = 0; is < 2; is++) {
			
				setup_FD3d<<<dimGrid,dimBlock>>>
					(SWFslice_d,v_d,ra_d,rb_d,rc_d,vel_d,ww, -caus*ab_h[4*is+0],ab_h[4*is+1],nx,ny);
				cusparseStatusX = cusparseCgtsvStridedBatch(cusparseHandleX,nx,rc_d,ra_d,rb_d,v_d,ny,nx);
				copy_transp_wfld3d_for<<<dimGrid,dimBlock>>>(v_d,SWFslice_d,nx,ny);
				//transposeDiagonal<<<dimGrid2,dimBlock2>>>(v_d,SWFslice_d,ny,nx);
				
				setup_FD3d<<<dimGrid,dimBlock>>>
					(SWFslice_d,v_d,ra_d,rb_d,rc_d,vel_d,ww, -caus*ab_h[4*is+2],ab_h[4*is+3],ny,nx);
				cusparseStatusY = cusparseCgtsvStridedBatch(cusparseHandleY,ny,rc_d,ra_d,rb_d,v_d,nx,ny);
				copy_transp_wfld3d_adj<<<dimGrid,dimBlock>>>(v_d,SWFslice_d,nx,ny);
				//transposeDiagonal<<<dimGrid2,dimBlock2>>>(v_d,SWFslice_d,nx,ny);

				setup_FD3d<<<dimGrid,dimBlock>>>
					(RWFslice_d,v_d,ra_d,rb_d,rc_d,vel_d,ww,-acaus*ab_h[4*is+0],ab_h[4*is+1],nx,ny);
				cusparseStatusX = cusparseCgtsvStridedBatch(cusparseHandleX,nx,rc_d,ra_d,rb_d,v_d,ny,nx);
				copy_transp_wfld3d_for<<<dimGrid,dimBlock>>>(v_d,RWFslice_d,nx,ny);
				//transposeDiagonal<<<dimGrid2,dimBlock2>>>(v_d,RWFslice_d,ny,nx);
				
				setup_FD3d<<<dimGrid,dimBlock>>>
					(RWFslice_d,v_d,ra_d,rb_d,rc_d,vel_d,ww,-acaus*ab_h[4*is+2],ab_h[4*is+3],ny,nx);
				cusparseStatusY = cusparseCgtsvStridedBatch(cusparseHandleY,ny,rc_d,ra_d,rb_d,v_d,nx,ny);
				copy_transp_wfld3d_adj<<<dimGrid,dimBlock>>>(v_d,RWFslice_d,nx,ny);
				//transposeDiagonal<<<dimGrid2,dimBlock2>>>(v_d,RWFslice_d,nx,ny);
				
			} /* END FD STEPS */

			/* High-angle filter FFT to kx */
			cufftExecC2C(plan, SWFslice_d, SWFslice_d, CUFFT_FORWARD);
			cufftExecC2C(plan, RWFslice_d, RWFslice_d, CUFFT_FORWARD);
	
			/* Fourier domain filtering */
			phase_correction3d<<<dimGrid, dimBlock>>>(SWFslice_d,minvel,ww,kxmap_d,kymap_d, caus,nx,ny,dz);
			phase_correction3d<<<dimGrid, dimBlock>>>(RWFslice_d,minvel,ww,kxmap_d,kymap_d,acaus,nx,ny,dz);
	
			/* High-angle filter FFT to kx */
			cufftExecC2C(plan, SWFslice_d, SWFslice_d, CUFFT_INVERSE);
			cufftExecC2C(plan, RWFslice_d, RWFslice_d, CUFFT_INVERSE);

			/* WAVEFIELD TAPER */
			wfld_taper3d<<<dimGrid, dimBlock>>>(SWFslice_d,tap_d,nx,ny);
			wfld_taper3d<<<dimGrid, dimBlock>>>(RWFslice_d,tap_d,nx,ny);
	
			/* XIG Imaging Condition */
 			imaging_condition<<<dimGrid,dimBlock>>>(SWFslice_d,RWFslice_d,xig_d,nx,ny,DNHX,DNHY);
		      
		    /* Insert slice */
			insert_slice<<<dimGrid, dimBlock>>>(swf_d,SWFslice_d,iw,nx,ny);
			insert_slice<<<dimGrid, dimBlock>>>(rwf_d,RWFslice_d,iw,nx,ny);
		      
		} /* END FREQUENCY LOOP */
				
		/* Copy back XIG slice and write output */
		cudaMemcpy( xig_h, xig_d, nx*ny*DNHX*DNHY*sizeof(float), cudaMemcpyDeviceToHost );
		
		sf_floatwrite(xig_h, nx*ny, Fxigo);
		
    } /* END DEPTH */

  if (wantwf) {
	/* Copy back SWF and RWF and write output */
    sf_warning("WRITING OUT RWF and SWF");
	cudaMemcpy( swf_h, swf_d, nx*ny*nw*sizeof(cuComplex), cudaMemcpyDeviceToHost );
	sf_complexwrite( swf_h, nx*nw, Fswfo);
	cudaMemcpy( rwf_h, rwf_d, nx*ny*nw*sizeof(cuComplex), cudaMemcpyDeviceToHost );
	sf_complexwrite( rwf_h, nx*ny*nw, Frwfo);
  }

  free(vel_h);   free(xig_h);
  free(swf_h);   free(rwf_h);
  free(kxmap_h); free(kymap_h);
  free(tap_h);   free(ab_h);
  
  cufftDestroy(plan);
  cusparseDestroy(cusparseHandleX);cusparseDestroy(cusparseHandleY);
  cudaFree(vel_d); cudaFree(xig_d); cudaFree(swf_d); cudaFree(rwf_d);
  cudaFree(tap_d); cudaFree(kxmap_d); cudaFree(kymap_d); 
  cudaFree(SWFslice_d); cudaFree(RWFslice_d);
  cudaFree(ra_d); cudaFree(rb_d); cudaFree(rc_d); cudaFree(v_d); 
  
  sf_close();
  exit(0);
}
