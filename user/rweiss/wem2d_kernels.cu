/* GPU Kernel Functions used in sfwem2d_gpu */

/*
  Authors: Jeffrey Shragge

  This file contains the GPU kernel functions called in the ewefd2d_gpu module from
  the Madagascar software package (http://www.reproducilitibly.org).  The calling
  functions for these kernels can be found in the file Mewefd2d_gpu.cu.  For more 
  information, see (Weiss and Shragge, "Solving 3D Anisotropic Elastic Wave 
  Equations on Parallel GPU Devices", GEOPHYSICS. http://software.seg.org/2012/0063)
*/


/*
  Copyright (C) 2012 University of Western Australia
  
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


#include <cuComplex.h>

#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define MAX(x, y) (((x) > (y)) ? (x) : (y))


// finite difference stencil coefficients are stored in constant device memory
__device__ __constant__ float C[4] = {0.040315157, 0.873981642, 0.457289566, 0.222691983};
__device__ __constant__ float trick = 0.11;

/*
 * 
 * GLOBAL KERNEL - wfld_taper:
 *   	Apply a wavefield taper to damp down boundary reflections
 *
 * launch configuration:
 *     	one thread per gridpoint
 * 		Grid	- (16,16,1)
 *		Block	- (nx/16,nw/16,1)
 *
 * Input:
 *		1) cuComplex *wf_d - An input wavefield array nx*nw
 *		2) float tap_d - Predefined taper array nx
 *		3) int nx - Number of x samples
 *		4) int nw - Number of frequencies
 * Output:
 * 
 * Return:
 * 		void
*/
__global__ void wfld_taper(	cuComplex *wf_d, 
				cuComplex *tap_d, 
				int nx, 
				int nw) 
/*< Taper wavefield >*/
{

	int ix = threadIdx.x + blockIdx.x * blockDim.x; /* X is ix */
	int iw = threadIdx.y + blockIdx.y * blockDim.y; /* Y is iw */

	if (ix < nx && iw < nw ) { /* If in Grid */
		int addr = iw*nx+ix;
		wf_d[addr] = cuCmulf( wf_d[addr] , tap_d[ix] );
	}
}


/*
 * 
 * GLOBAL KERNEL - prop_split_step:
 *   	Propagation for zero wavenumber
 *
 * launch configuration:
 *     	one thread per gridpoint
 * 		Grid	- (16,16,1)
 *		Block	- (nx/16,nw/16,1)
 *
 * Input:
 *
 * Output:
 * 
 * Return:
 * 		void
*/
__global__ void prop_SSF(	cuComplex *wld_d, 
							float *vel_d, 
							float dw,
							float ow,
							float weisign, 
							int iz,
							float dz,
							int nx,
							int nw) 
/*< Split-step bulk phase propagator >*/							
{
	int ix = threadIdx.x + blockIdx.x * blockDim.x; /* X is ix */
	int iw = threadIdx.y + blockIdx.y * blockDim.y; /* Y is iw */
	float ww = 2.f*SF_PI*((float)iw*dw+ow);
	
	if (ix < nx && iw < nw ) {
		int vaddr = iz*nx+ix; /* Velocity address */
		int  addr = iw*nx+ix; /* wavefield address */
		float tmp = - weisign * dz * ww / vel_d[vaddr];
		cuComplex phsft = make_cuFloatComplex(cos(tmp),sin(tmp));
		wld_d[addr] = cuCmulf(wld_d[addr],phsft);
	}
}


/*
 * 
 * KERNEL - prop_finite_difference:
 *   	Propagation for zero wavenumber
 *
 * launch configuration:
 *     	one thread per gridpoint
 * 		Grid	-
 *		Block	-
 *
 * Input:
 *
 * Output:
 * 
 * Return:
 * 		void
*/
__global__ void setup_FD(cuComplex *wld_d,
			 cuComplex *v_d,
			 cuComplex *ra_d,
			 cuComplex *rb_d,
			 cuComplex *rc_d, 
			 float *vel_d, 
			 float dw,
			 float ow,
			 float ax,
			 float bx, 
			 int nx, 
			 int nw,
			 int iz)
/*< Finite-difference propagator >*/						
{
  
  int ix = threadIdx.x + blockIdx.x * blockDim.x; /* X is ix */
  int iw = threadIdx.y + blockIdx.y * blockDim.y; /* Y is iw */
  
  cuComplex url,urr,ror,rol;
  cuComplex zero = make_cuFloatComplex(0.f,0.f);
  cuComplex one  = make_cuFloatComplex(1.f,0.f);
  cuComplex two  = make_cuFloatComplex(2.f,0.f);
  
  if (iw < nw && ix < nx) { /* If in range */
    int vaddr = iz*nx+ix;
    int  addr = iw*nx+ix;
    float ca =  vel_d[vaddr] / (2.f*SF_PI*((float)iw*dw+ow));
    ror = make_cuFloatComplex(trick + bx*ca*ca , ax*ca);
    rol = cuConjf(ror);		
    
    if (ix == 0) { /* CASE WHERE IX=0  */
      url = zero; 
      urr = wld_d[addr+1];
      v_d [addr] = cuCaddf(cuCmulf(wld_d[addr],cuCsubf(one,cuCmulf(two,ror))),cuCmulf(ror,(cuCaddf(url,urr))));
      ra_d[addr] = cuCsubf(one , cuCmulf(two,rol));
      rb_d[addr] = rol;
      rc_d[addr] = zero;
    } else if (ix == nx-1)  { /* CASE WHERE IX=NX-1  */
      url = wld_d[addr-1];
      urr = zero;
      v_d [addr] = cuCaddf(cuCmulf(wld_d[addr],cuCsubf(one,cuCmulf(two,ror))),cuCmulf(ror,(cuCaddf(url,urr))));
      ra_d[addr] = cuCsubf(one , cuCmulf(two,rol));
      rb_d[addr] = zero;
      rc_d[addr] = rol;
    } else { /* NON-END CASE */
      url = wld_d[addr-1];
      urr = wld_d[addr+1];
      v_d [addr] = cuCaddf(cuCmulf(wld_d[addr],cuCsubf(one,cuCmulf(two,ror))),cuCmulf(ror,(cuCaddf(url,urr))));
      ra_d[addr] = cuCsubf(one , cuCmulf(two,rol));
      rb_d[addr] = rol;
      rc_d[addr] = rol;			
    }
    
  } /* END in range */
  
}


/*
 * 
 * KERNEL - phase_correction:
 *   	Propagation for zero wavenumber
 *
 * launch configuration:
 *     	one thread per gridpoint
 * 		Grid	-
 *		Block	-
 *
 * Input:
 *
 * Output:
 * 
 * Return:
 * 		void
 */
__global__ void phase_correction( 	cuComplex *wfl_d,
					float vmin,
					float dw,
					float ow,
					float *kxmap_d,
					int iz,
					float weisign,
					int nkx,
					int nw,
					float dz) 
/*< Phase shift correction >*/
{
  int ikx = threadIdx.x + blockIdx.x * blockDim.x; /* X is ix */
  int iw  = threadIdx.y + blockIdx.y * blockDim.y; /* Y is iw */
  
  /* Ensure in the correct location */
  if (ikx < nkx && iw < nw ) {	
    cuComplex zero = make_cuFloatComplex(0.f,0.f);						
    float w_v0 = 2.f*SF_PI*((float)iw*dw+ow) / vmin;
    float avX2 = 1.f/(w_v0*w_v0);
    int attr = iw*nkx+ikx;
    float kx2 = pow(kxmap_d[ikx],2);
    float sr2 = kx2 * avX2;	
    
    /* Filter out evanescent energy */
    if (sr2 > 1.f ) {
      wfl_d[attr]=zero;
    } else { /* Apply filter */
      float tx = kx2*avX2;
      float phsft=-dz*w_v0*( sqrt(1.f-sr2)-(1.f-C[0]*tx/(1.f-C[1]*tx)-C[2]*tx/(1.f-C[3]*tx)));
      wfl_d[attr] = cuCmulf(wfl_d[attr],make_cuFloatComplex(cos(phsft),sin(phsft)*weisign));
    }
    
  }								
  
}

__global__ void copy_wfld( 	cuComplex *in, 
				cuComplex *out, 
				int nx, 
				int nw)
/*< Copy over a wavefield from one location to the next >*/
{
  int ix = threadIdx.x + blockIdx.x * blockDim.x; /* X is ix */
  int iw = threadIdx.y + blockIdx.y * blockDim.y; /* Y is iw */
  
  if (ix < nx && iw < nw ) {
    int attr = iw * nx + ix;
    out[attr] = in[attr];
  }	
  
}

__global__ void rescale_wfld(cuComplex *wld, int nx, int nw)
/*< Copy over a wavefield from one location to the next >*/
{
  int ix = threadIdx.x + blockIdx.x * blockDim.x; /* X is ix */
  int iw = threadIdx.y + blockIdx.y * blockDim.y; /* Y is iw */
  
  if (ix < nx && iw < nw ) {
    int attr = iw * nx + ix;
    wld[attr] = cuCdivf(wld[attr],make_cuFloatComplex((float)nx ,0.f) ) ;
  }	
  
}

/*
 * 
 * KERNEL - illumination:
 *   	Compute the auto-correction of the source wavefield
 *
 * launch configuration:
 *     	one thread per gridpoint
 * 		Grid	-
 *		Block	-
 *
 * Input:
 *		wf
 * Output:
 * 
 * Return:
 * 		void
*/
__global__ void illumination(cuComplex *wf_d, 
				  cuComplex *illum_d, 
				  int iz,
				  int nx, 
				  int nw)
/*< Illumination (diagonal) */
{
  int ix = threadIdx.x + blockIdx.x * blockDim.x; /* X is ix */
  
  if ( ix < nx ) {
    
    cuComplex sum=make_cuFloatComplex(0.f,0.f);
    for (int iw=0; iw < nw; iw++) {
      int waddr = iw*nx+ix;
      sum = cuCaddf(sum, cuCmulf(cuConjf(wf_d[waddr]),wf_d[waddr]));
    }
    illum_d[ix]=sum;
   
  } // END NX LOOP
}


/*
 * 
 * KERNEL - imaging_condition:
 *   	Cross correlation of SWF and RWF (no shift)
 *
 * launch configuration:
 *     	one thread per gridpoint
 * 		Grid	-
 *		Block	-
 *
 * Input:
 *
 * Output:
 * 
 * Return:
 * 		void
 */
__global__ void imaging_condition(cuComplex *swf_d, 
				  cuComplex *rwf_d, 
				  cuComplex *img_d, 
				  int iz,
				  int nh,
				  int nx, 
				  int nw)
/*< Imaging Condition with Parallel reduction */
{
  int ix = threadIdx.x + blockIdx.x * blockDim.x; /* X is ix */
  
  if ( ix < nx ) {
    
    cuComplex sum=make_cuFloatComplex(0.f,0.f);
    for (int iw=0; iw < nw; iw++) {
      int waddr = iw*nx+ix;
      sum = cuCaddf(sum, cuCmulf(cuConjf(swf_d[waddr]),rwf_d[waddr]));
    }
    img_d[ix]=sum;
   
  } // END NX LOOP
}


/*
 * 
 * KERNEL - xig_condition:
 *   	Cross correlation of SWF and RWF (with horizontal shift)
 *
 * launch configuration:
 *     	one thread per gridpoint
 * 		Grid	-
 *		Block	-
 *
 * Input:
 *
 * Output:
 * 
 * Return:
 * 		void
 */
__global__ void xig_condition(cuComplex *swf_d, 
			      cuComplex *rwf_d, 
			      cuComplex *img_d, 
			      int iz,
			      int nh,
			      int nx, 
			      int nw)
/*< XIG Imaging Condition with Shifted correlation */
{
  int ix = threadIdx.x + blockIdx.x * blockDim.x; /* X is ix */
  int iaddr,waddr_N,waddr_P;

  if ( ix < nx ) {
    
    for (int ih = -(int)(nh-1)/2; ih < (int)(nh+1)/2; ih++) {
    	cuComplex sum=make_cuFloatComplex(0.f,0.f);
      
      if ( ix+ih > -1 && ix-ih > -1 && ix-ih < nx && ix+ih < nx) {
		waddr_N= ix-ih-nx;
		waddr_P= ix+ih-nx;
		for (int iw=0; iw < nw; iw++) {
	  		waddr_N += nx;
	  		waddr_P += nx;
		    sum = cuCaddf(sum,cuCmulf(cuConjf(swf_d[waddr_N]),rwf_d[waddr_P]));
		} // END NW LOOP
		iaddr = ix*nh+ih+(int)(nh-1)/2;
		img_d[iaddr]=sum;	
      } // END IF LOOP
      
    } // END NH LOOP
    
  } // END NX LOOP
  
}

/*
 * 
 * KERNEL - make_adj_wflds:
 *   	Make adjoint SWF and RWF (no coupling)
 *      SWFadj = conj(image)*RWF
 *      RWFadj =      image *SWF
 *
 * launch configuration:
 *     	one thread per gridpoint
 * 		Grid	-
 *		Block	-
 *
 * Input:
 *
 * Output:
 * 
 * Return:
 * 		void
 */
__global__ void make_adj_wflds( cuComplex *xig_d,
				cuComplex *swf_d,
				cuComplex *rwf_d,
				cuComplex *swfadj_d,
				cuComplex *rwfadj_d,
				int nh,
				int nx,
				int nw)
/*< Make the adjoint wavefields >*/
{
  
  int ix = threadIdx.x + blockIdx.x * blockDim.x; /* X is ix */
  int iw = threadIdx.y + blockIdx.y * blockDim.y; /* Y is iw */
  if (ix < nx && iw < nw) { /* IN GRID */
    
    for (int ih= -(int)(nh-1)/2; ih < (int)(nh+1)/2; ih++) { /* SWFADJ */
      
      if ( ix+ih > 0 && ix-ih > 0 && ix-ih < nx && ix+ih < nx) {
		int iaddr = ix*nh+ih+(int)(nh-1)/2;
		int mind = iw*nx + ix-ih;
		int pind = iw*nx + ix+ih;
		cuComplex xx=xig_d[iaddr];
		swfadj_d[pind] = cuCaddf(swfadj_d[pind],cuCmulf(cuConjf(xx),rwf_d[mind]));
		rwfadj_d[mind] = cuCaddf(rwfadj_d[mind],cuCmulf(        xx ,swf_d[pind]));
      } 
    }
    
  }
  
}

/*
 * 
 * KERNEL - make_adj_wflds_coupled:
 *   	Make adjoint wavefields [as1,ar1] and [as2,ar2] (coupling)
 *      penDSO = |h|
 *      pen1 = sech(r1^2)
 *      pen2 = sech(r2^2)
 *      as1 = conj(r1)*ur1*(epsDSO*penDSO^2+eps4D*P2^2+eps4D*r2^2*tanh(r2^2)*P2^2)
 *      ar1 =      r1 *us1*(epsDSO*penDSO^2+eps4D*P2^2+eps4D*r2^2*tanh(r2^2)*P2^2)
 *      as2 = conj(r2)*ur2*(epsDSO*penDSO^2+eps4D*P2^2+eps4D*r2^2*tanh(r2^2)*P2^2)
 *      ar2 =      r2 *us2*(epsDSO*penDSO^2+eps4D*P2^2+eps4D*r2^2*tanh(r2^2)*P2^2) 
 *
 * launch configuration:
 *     	one thread per gridpoint
 * 		Grid	-
 *		Block	-
 *
 * Input:
 *
 * Output:
 * 
 * Return:
 * 		void
 */
__global__ void make_adj_wflds_coupled( cuComplex *xig1_d,
				cuComplex *xig2_d,
				cuComplex *us1_d,
				cuComplex *ur1_d,
				cuComplex *as1_d,
				cuComplex *ar1_d,
				cuComplex *us2_d,
				cuComplex *ur2_d,
				cuComplex *as2_d,
				cuComplex *ar2_d,
				float epsDSO,
				int hzero,
				float eps4D,
				float epsNORM,
				float dx,
				int nh,
				int nx,
				int nw)
/*< Make the adjoint wavefields >*/
{
  
  int ix = threadIdx.x + blockIdx.x * blockDim.x; /* X is ix */
  int iw = threadIdx.y + blockIdx.y * blockDim.y; /* Y is iw */
  if (ix < nx && iw < nw) { /* IN GRID */
    
    for (int ih= -(int)(nh-1)/2; ih < (int)(nh+1)/2; ih++) { /* H-LOOP */
  
      if ( ix+ih > 0 && ix-ih > 0 && ix-ih < nx && ix+ih < nx) {
		int iaddr = ix*nh+ih+(int)(nh-1)/2;		
		/* int  ind  = iw*nx+ix   ; */
		int mind  = iw*nx+ix-ih;
		int pind  = iw*nx+ix+ih;
		
		/* Local image variables */
		register cuComplex x1 = xig1_d[iaddr];
		register cuComplex x2 = xig2_d[iaddr];
		
		cuComplex us11 = us1_d[mind];
		cuComplex us22 = us2_d[mind];

		cuComplex ur11 = ur1_d[pind];		
		cuComplex ur22 = ur2_d[pind];
		
	 	x1  = cuCmulf(x1,cuCmulf(cuConjf(us11),ur11));
		x2  = cuCmulf(x2,cuCmulf(cuConjf(us22),ur22));
		
		/* Baseline adjoint source terms */
		as1_d[mind] = cuCaddf(as1_d[mind],cuCmulf(cuConjf(x1),ur11));
  		ar1_d[pind] = cuCaddf(ar1_d[pind],cuCmulf(        x1 ,us11));
				
		/* Monitor adjoint source terms */
		as2_d[mind] = cuCaddf(as2_d[mind],cuCmulf(cuConjf(x2),ur22));
 		ar2_d[pind] = cuCaddf(ar2_d[pind],cuCmulf(        x2 ,us22));		
		
      }  /* END IF */
    } /* END OFFSET ih */

  } /* END IF IN GRID */
}

/*
 * 
 * KERNEL - compute:
 *   	Compute gradient by cross correlating 
 * 	      wavefields and adjoint wavefields
 *
 *      gradient = conj(SWFadj)*SWF+conj(RWFadj)*RWF
 *
 * launch configuration:
 *     	one thread per gridpoint
 * 		Grid	-
 *		Block	-
 *
 * Input:
 *
 * Output:
 * 
 * Return:
 * 		void
 */
__global__ void compute_gradient(float *grd_d,
				 cuComplex *us_d,
				 cuComplex *ur_d,
				 cuComplex *as_d,
				 cuComplex *ar_d,
				 int nx,
				 float ow,
				 float dw,
				 int nw)
/*< Compute the gradient >*/
{
  int ix = threadIdx.x + blockIdx.x * blockDim.x; /* X is ix */
  if (ix < nx ) { /* IF IN GRID */
    float sum = 0.f;
    for (int iw=0; iw<nw; iw++) {
      int wa = iw*nx+ix;
      //float ww = (float)(iw-1.f)*dw+ow;
      //sum+=cuCimagf(cuCmulf(ur_d[wa],cuConjf(ar_d[wa])));// - cuCimagf(cuCmulf(ar_d[wa],cuConjf(ur_d[wa])));
      sum+=cuCimagf(cuCaddf(cuCmulf((us_d[wa]),cuConjf(as_d[wa])),cuCmulf((ur_d[wa]),cuConjf(ar_d[wa]))));
    }
    grd_d[ix] = sum;
  }
}

