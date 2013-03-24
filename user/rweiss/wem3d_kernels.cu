/* GPU Kernel Functions used in sfwem3d_gpu */

/*
  Authors: Jeffrey Shragge

  This file contains the GPU kernel functions called in the ewefd2d_gpu module from
  the Madagascar software package (http://www.reproducilitibly.org).  The calling
  functions for these kernels can be found in the file Mewefd2d_gpu.cu.  For more 
  information, see (Weiss and Shragge, "Solving 3D Anisotropic Elastic Wave 
  Equations on Parallel GPU Devices", GEOPHYSICS. http://software.seg.org/2012/0063)
*/


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

#include <cuComplex.h>

#define DNHX 1
#define DNHY 1
#define BLOCKX 16
#define BLOCKY 16
// finite difference stencil coefficients are stored in constant device memory
__device__ __constant__ float C[4] = {0.040315157, 0.873981642, 0.457289566, 0.222691983};
__device__ __constant__ float trick = 0.11;

/*
 * 
 * GLOBAL KERNEL - wfld_taper3d:
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
 *		4) int ny - Number of y samples
 * Output:
 * 
 * Return:
 * 		void
*/
__global__ void wfld_taper3d(	
 				cuComplex *wf_d, 
				cuComplex *tap_d, 
				int nx, 
				int ny) 
/*< Taper wavefield >*/
{
	unsigned int ix = threadIdx.x + blockIdx.x * blockDim.x; /* X is ix */
	unsigned int iy = threadIdx.y + blockIdx.y * blockDim.y; /* Y is iy */

	if (ix < nx && iy < ny ) { /* If in Grid */
		int addr = iy*nx+ix;
		wf_d[addr] = cuCmulf( wf_d[addr] , tap_d[addr] );
	}
}

/*
 * 
 * GLOBAL KERNEL - prop_SSF3d:
 *   	Propagation for zero wavenumber
 *
 * launch configuration:
 *     	one thread per gridpoint
 * 		Grid	- (16,16,1)
 *		Block	- (nx/16,ny/16,1)
 *
 * Input:
 *
 * Output:
 * 
 * Return:
 * 		void
*/
__global__ void prop_SSF3d(	cuComplex *wld_d, 
							float *vel_d, 
							float ww,
							float weisign, 
							float dz,
							int nx,
							int ny) 
/*< Split-step bulk phase propagator >*/							
{
	unsigned int ix = threadIdx.x + blockIdx.x * blockDim.x; /* X is ix */
	unsigned int iy = threadIdx.y + blockIdx.y * blockDim.y; /* Y is iy */
	
	if (ix < nx && iy < ny ) { /* If in grid */
		unsigned int addr = iy*nx+ix; /* wavefield address */
		float tmp = - weisign * dz * ww / vel_d[addr];
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
__global__ void setup_FD3d(cuComplex *wld_d,
						 cuComplex *v_d,
						 cuComplex *ra_d,
						 cuComplex *rb_d,
						 cuComplex *rc_d, 
						 float *vel_d, 
						 float ww,
						 float ax,
						 float bx, 
						 int nx, 
						 int ny)
/*< Finite-difference propagator >*/						
{
	unsigned int ix = threadIdx.x + blockIdx.x * blockDim.x; /* X is ix */
	unsigned int iy = threadIdx.y + blockIdx.y * blockDim.y; /* Y is iy */

	cuComplex url,urr,ror,rol;
	cuComplex zero = make_cuFloatComplex(0.f,0.f);
	cuComplex one  = make_cuFloatComplex(1.f,0.f);
	cuComplex two  = make_cuFloatComplex(2.f,0.f);
	
	if (iy < ny && ix < nx) { /* If in range */
		unsigned int addr = iy*nx+ix;
		float ca = vel_d[addr] / ww;
		ror = make_cuFloatComplex(trick + bx*ca*ca , ax*ca);
		rol = cuConjf(ror);		

		if (ix == 0) { /* CASE WHERE IX=0  */
			url = zero; 
			urr = wld_d[addr+1];
			v_d [addr] = cuCaddf(cuCmulf(wld_d[addr],cuCsubf(one,cuCmulf(two,ror))),cuCmulf(ror,(cuCaddf(url,urr))));
			ra_d[addr] = cuCsubf(one,cuCmulf(two,rol));
			rb_d[addr] = rol;
			rc_d[addr] = zero;
			
		} else if (ix == nx-1)  { /* CASE WHERE IX=NX-1  */
			url = wld_d[addr-1];
			urr = zero;
			v_d [addr] = cuCaddf(cuCmulf(wld_d[addr],cuCsubf(one,cuCmulf(two,ror))),cuCmulf(ror,(cuCaddf(url,urr))));
			ra_d[addr] = cuCsubf(one,cuCmulf(two,rol));
			rb_d[addr] = zero;
			rc_d[addr] = rol;
			
		} else { /* NON-END CASE */
			url = wld_d[addr-1];
			urr = wld_d[addr+1];
			v_d [addr] = cuCaddf(cuCmulf(wld_d[addr],cuCsubf(one,cuCmulf(two,ror))),cuCmulf(ror,(cuCaddf(url,urr))));
			ra_d[addr] = cuCsubf(one,cuCmulf(two,rol));
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
__global__ void phase_correction3d(cuComplex *wfl_d,
									float vmin,
									float ww,
									float *kxmap_d,
									float *kymap_d,
									float weisign,
									int nkx,
									int nky,
									float dz) 
/*< Phase shift correction >*/
{
	unsigned int ikx = threadIdx.x + blockIdx.x * blockDim.x; /* X is ix */
	unsigned int iky = threadIdx.y + blockIdx.y * blockDim.y; /* Y is iy */
  	cuComplex zero = make_cuFloatComplex(0.f,0.f);						
	cuComplex dimnorm = make_cuFloatComplex( 1.f/(float)(nkx*nky),0.f);
	
  	/* Ensure in the correct location */
  	if (ikx < nkx && iky < nky) {	
		float w_v0 = ww / vmin;
		float avX2 = 1.f/(w_v0*w_v0);
		float kx2 = pow(kxmap_d[ikx],2);
		float ky2 = pow(kymap_d[iky],2);
	    float sr2 = (kx2 + ky2) * avX2;	
	        	
    	/* Filter out evanescent energy */
    	if (sr2 > 1.f ) {
	    	unsigned int attr = iky*nkx+ikx;
    	  	wfl_d[attr]=zero;
    	} else { /* Apply filter */
        	unsigned int attr = iky*nkx+ikx;
		  	float tx = kx2*avX2;
    	  	float ty = ky2*avX2;
    	  	float t1 = -C[0]*tx/(1.f-C[1]*tx);
    	  	float t2 = -C[2]*tx/(1.f-C[3]*tx);
    	  	float t3 = -C[0]*ty/(1.f-C[1]*ty);
			float t4 = -C[2]*ty/(1.f-C[3]*ty);
    	  	float phsft=-dz*w_v0*( sqrt(1.f-sr2)-(1.f+t1+t2+t3+t4));
    	  	wfl_d[attr] = cuCmulf(cuCmulf(wfl_d[attr],dimnorm),make_cuFloatComplex(cos(phsft),sin(phsft)*weisign));
  		}							
	}															
}

__global__ void extract_slice(cuComplex *volume,
							  cuComplex *slice,
							  int iw,
							  int nx,
							  int ny)
/*< Extract 2D slice from 3D volume >*/
{
	unsigned int ix = threadIdx.x + blockIdx.x * blockDim.x; /* X is ix */
	unsigned int iy = threadIdx.y + blockIdx.y * blockDim.y; /* Y is iy */
	if (ix < nx && iy < ny ) {
		unsigned int sattr =    iy*nx + ix;
		unsigned int vattr = iw*ny*nx + sattr;
		slice[sattr] = volume[vattr];
	}
}							  

__global__ void insert_slice(cuComplex *volume,
							 cuComplex *slice,
							  int iw,
							  int nx,
							  int ny)
/*< Insert 2D slice into 3D volume >*/
{
	unsigned int ix = threadIdx.x + blockIdx.x * blockDim.x; /* X is ix */
	unsigned int iy = threadIdx.y + blockIdx.y * blockDim.y; /* Y is iy */
	if (ix < nx && iy < ny ) {
		unsigned int sattr =    iy*nx + ix;
		unsigned int vattr = iw*ny*nx + sattr;
		volume[vattr] = slice[sattr];
	}
}							  


__global__ void transposeDiagonal(cuComplex *odata,  
            cuComplex *idata, int width, int height) 
{ 
  __shared__ cuComplex tile[32][32+1]; 
  int blockIdx_x, blockIdx_y; 
  // diagonal reordering 
  if (width == height) { 
    blockIdx_y = blockIdx.x; 
    blockIdx_x = (blockIdx.x+blockIdx.y)%gridDim.x;
  } else { 
    int bid = blockIdx.x + gridDim.x*blockIdx.y; 
    blockIdx_y = bid%gridDim.y; 
    blockIdx_x = ((bid/gridDim.y)+blockIdx_y)%gridDim.x; 
  }     
  int xIndex = blockIdx_x*32 + threadIdx.x; 
  int yIndex = blockIdx_y*32 + threadIdx.y;  
  int index_in = xIndex + (yIndex)*width; 
  xIndex = blockIdx_y*32 + threadIdx.x; 
  yIndex = blockIdx_x*32 + threadIdx.y; 
  int index_out = xIndex + (yIndex)*height; 
  for (int i=0; i<32; i+=8) { 
    tile[threadIdx.y+i][threadIdx.x] = idata[index_in+i*width]; 
  } 
   
  __syncthreads(); 
  for (int i=0; i<32; i+=8) { 
    odata[index_out+i*height] = tile[threadIdx.x][threadIdx.y+i]; 
  } 

} 

__global__ void copy_transp_wfld3d_for(cuComplex *idata, cuComplex *odata, int width, int height)
{
	__shared__ cuComplex block[BLOCKX][BLOCKX+1];
	
	// read the matrix tile into shared memory
	unsigned int xIndex = blockIdx.x * BLOCKX + threadIdx.x;
	unsigned int yIndex = blockIdx.y * BLOCKX + threadIdx.y;
	if((xIndex < width) && (yIndex < height))
	{
		unsigned int index_in = yIndex * width + xIndex;
		block[threadIdx.y][threadIdx.x] = idata[index_in];
	}

	__syncthreads();
	// write the transposed matrix tile to global memory
	xIndex = blockIdx.y * BLOCKX + threadIdx.x;
	yIndex = blockIdx.x * BLOCKX + threadIdx.y;
	if((xIndex < height) && (yIndex < width))
	{
		unsigned int index_out = yIndex * height + xIndex;
		odata[index_out] = block[threadIdx.x][threadIdx.y];
	}
}

__global__ void copy_transp_wfld3d_adj(cuComplex *idata, cuComplex *odata, int width, int height)
{
	__shared__ cuComplex block[BLOCKX][BLOCKX+1];
	
	// read the matrix tile into shared memory
	unsigned int xIndex = blockIdx.y * BLOCKX + threadIdx.x; 
	unsigned int yIndex = blockIdx.x * BLOCKX + threadIdx.y; 
	
	if((xIndex < height) && (yIndex < width))
	{
		unsigned int index_in = yIndex * height + xIndex; 
		block[threadIdx.y][threadIdx.x] = idata[index_in];
	}

	__syncthreads();
	// write the transposed matrix tile to global memory
	xIndex = blockIdx.x * BLOCKX + threadIdx.x;
	yIndex = blockIdx.y * BLOCKX + threadIdx.y;
	if((xIndex < width) && (yIndex < height))
	{
		unsigned int index_out =  yIndex * width + xIndex;
		odata[index_out] = block[threadIdx.x][threadIdx.y];
	}
}



__global__ void imaging_condition(cuComplex *swf_d, 
								  cuComplex *rwf_d, 
								  float *img_d, 
								  int nx,
								  int ny, 
								  int nhx,
								  int nhy)
/*< Imaging Condition */
{
  	cuComplex zero = make_cuFloatComplex(0.f,0.f);						
	unsigned int ix = threadIdx.x + blockIdx.x * blockDim.x; /* X is ix */
	unsigned int iy = threadIdx.y + blockIdx.y * blockDim.y; /* Y is iy */
	
	if ( ix < nx && iy < ny ) {

		if (DNHX == 1 && DNHY == 1 ) {

			unsigned int addr = iy*nx+ix;
			img_d[addr] += cuCrealf(cuCmulf(cuConjf(swf_d[addr]),rwf_d[addr]));

		} else {
			int nsx = BLOCKX+2*DNHX;
			int nsy = BLOCKY+2*DNHY;
			__shared__ cuComplex lswf_d[(BLOCKX+2*DNHX)*(BLOCKY+2*DNHY)];
			__shared__ cuComplex lrwf_d[(BLOCKX+2*DNHX)*(BLOCKY+2*DNHY)];

			/* Bring wavefields into shared memory */
			/* Include NHX and NHY boundary */
			for (int isy=0; isy < nsy; isy++) {
				for (int isx=0; isx < nsx; isx++) {
					
					/* Local Address */
					unsigned int laddr = isy*nsx + isx;
					lswf_d[laddr]=zero;
					lrwf_d[laddr]=zero;
					
					if ( iy-nhy+isy > 0  && iy-nhy+isy < ny) {
						if ( ix-nhx+isx > 0 && ix-nhx+isx < nx) {
						
							/* Global Address */
							unsigned int gaddr = (iy-nhy) * nx + ix-nhx;
							
							/* Updated local wavefields */
							lswf_d[laddr] = swf_d[gaddr];
							lrwf_d[laddr] = rwf_d[gaddr];
							
						} /* End X IF */
					} /* End Y IF */
				} /* END X LOOP */
			} /* END Y LOOP */
		 
		 
		 	/* Imaging Condition */
	  		/* IHY LOOP */
			for (int ihy = -(int)((nhy-1)/2.f); ihy < (int)((nhy-1)/2.f); ihy++) {	
				/* IHX LOOP */
				for (int ihx = -(int)((nhx-1)/2.f); ihx < (int)((nhx-1)/2.f); ihx++) {
					unsigned int waddr_N = (iy-ihy)*nsx + ix-ihx;
					unsigned int waddr_P = (iy+ihy)*nsx + ix+ihx;
					
					/* Global image output nx-ny-nhx-nhy-nz */
					unsigned int iaddr = ihy*nhx*ny*nx+ihx*ny*nx+iy*nx+ix;
					
					img_d[iaddr] += cuCrealf(cuCmulf(cuConjf(lswf_d[waddr_N]),lrwf_d[waddr_P]));
	
				} /* END NHX LOOP */
			} /* END NHY LOOP */
		} /* END Zero OFFSET LOOP */			
	} /* END CHECK LOOP */
		
}

__global__ void make_adj_wflds( float *xig_d,
								cuComplex *swf_d,
								cuComplex *rwf_d,
								cuComplex *swfadj_d,
								cuComplex *rwfadj_d,
								int nh,
								int nx,
								int nw)
/*< Make the adjoint wavefields >*/
{

	unsigned int ix = threadIdx.x + blockIdx.x * blockDim.x; /* X is ix */
	unsigned int iw = threadIdx.y + blockIdx.y * blockDim.y; /* Y is iw */
 	if (ix < nx && iw < nw) { /* IN GRID */
	  	
		for (int ih= -(int)(nh-1)/2.f; ih < (int)(nh-1)/2.f; ih++) { /* SWFADJ */
		
			if ( ix+ih > 0 && ix-ih > 0 && ix-ih < nx && ix+ih < nx) {
				unsigned int iaddr = ix*nh+ih+(int)(nh-1)/2.f;
	    		unsigned int mind = iw*nx + ix-ih;
	    		unsigned int pind = iw*nx + ix+ih;
	    		swfadj_d[pind] = cuCaddf(swfadj_d[pind],cuCmulf(make_cuFloatComplex(xig_d[iaddr],0.f),rwf_d[mind]));
			} 
		}
	
	  	for (int ih=-(int)(nh-1)/2.f; ih < (int)(nh-1)/2.f; ih++) { /* RWFADJ */
	  	
			if ( ix+ih > 0 && ix-ih > 0 && ix-ih < nx && ix+ih < nx) {
				unsigned int iaddr = ix*nh+ih+(int)(nh-1)/2.f;
	    		unsigned int mind = iw*nx + ix-ih;
	    		unsigned int pind = iw*nx + ix+ih;
		    	rwfadj_d[mind] = cuCaddf(rwfadj_d[mind],cuCmulf(make_cuFloatComplex(xig_d[iaddr],0.f),swf_d[pind]));
			} 
	 	}
	}

}

__global__ void compute_gradient(	float *grd_d,
									cuComplex *s_d,
									cuComplex *r_d,
									cuComplex *sa_d,
									cuComplex *ra_d,
									int nx,
									float ow,
									float dw,
									int nw)
/*< Compute the gradient >*/
{
	unsigned int ix = threadIdx.x + blockIdx.x * blockDim.x; /* X is ix */
 	if (ix < nx ) { /* IF IN GRID */
		float sum = 0.f;
		for (int iw=0; iw<nw; iw++) {
			unsigned int wa = iw*nx+ix;
			sum+=-cuCimagf(cuCaddf(cuCmulf(sa_d[wa],cuConjf(s_d[wa])),cuCmulf(r_d[wa],cuConjf(ra_d[wa]))));
		}
		grd_d[ix] = sum;
	}
}

