/* Wave-equation migration caller (single shot) routine for 2D ISOTROPIC FD migration */
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
#include <rsf.h>
#include <math.h>
#include <complex.h>
#ifdef _OPENMP
#include <omp.h>
#endif
/*^*/
#include "wem2d_iso_kernel.h"
#include "fft1_axis2.h"
#include "wem2d_iso_wemig.h"
/*^*/

static void wem2d_iso_rwetaper_init(void);
static void wem2d_iso_rwetaper(sf_complex **dax, int id);

static int nx; /* num inline */
static int nw; /* num frequencies */
static int nz; /* num depths */
static int nh; /* num xig shifts */
static int nxtap; /* taper length */

static float dx; /* inline sampling */
static float dz; /* depth sampling */
static float dw; /* frequency sampling */
static float dh; /* xig shift sampling */

static float ox; /* inline offset */
static float oz; /* depth offset */
static float ow; /* frequency offset */
static float oh; /* xig shift offset */

static float ****txig; /* Large temp xig volume */
static float *vmin; /* Min velocity as function of depth */
static float *tap; /* Boundary taper */
static sf_complex **sax; /* Single w   source wavefield (OMP) */
static sf_complex **rax; /* Single w receiver wavefield (OMP) */
static float *ww;

/* Fourier components */
static int nk; 
static float ok,dk;

static int nth;
/*-------------------------------------------------------------*/
/* DEFINES */

#define XLOOP(a) for (ix=0; ix<nx; ix++) { {a} }
#define ZLOOP(a) for (iz=0; iz<nz; iz++) { {a} }
#define ZXLOOP(a) for (iz = 0; iz<nz; iz++){ for (ix=0; ix<nx; ix++) { {a} }}
#define ZNXLOOP(a) for (iz = 0; iz<nz; iz++) {	\
	for (ith = 0; ith < nth; ith++) {	\
	    for (ix=0; ix<nx; ix++) {		\
		{a} }}}
#define ZNXHLOOP(a) for (iz = 0; iz<nz; iz++) {	\
	for (ith = 0; ith < nth; ith++) {	\
	    for (ix=0; ix<nx; ix++) {		\
		for (ih=0; ih<nh; ih++) {	\
		    {a} }}}}

/*-------------------------------------------------------------*/
/* Initialization routine */
void wem2d_iso_wemig_init(int  nx_in,  int nz_in,   int nw_in, int nh_in,
			  float ox_in, float oz_in, float ow_in, float oh_in,
			  float dx_in, float dz_in, float dw_in, float dh_in,
			  int nxtap_in, float **vel)
/*< initialize >*/
{
    int iz,ix;

    nx=nx_in; ox=ox_in; dx=dx_in;
    nz=nz_in; oz=oz_in; dz=dz_in; 
    nh=nh_in; oh=oh_in; dh=dh_in;
    nw=nw_in; ow=2.f*SF_PI*ow_in; dw=2.f*SF_PI*dw_in;
  
    nxtap=nxtap_in;
    nth = 1;

#ifdef _OPENMP
#pragma omp parallel
    {
	nth = omp_get_num_threads();
    }
#endif
    sf_warning("USING the following number of threads: %d",nth);
  
    /* ALLOCATIONS */
    txig = sf_floatalloc4(nh,nx,nth,nz);
    vmin = sf_floatalloc(nz);
    ww   = sf_floatalloc(nth);
    tap  = sf_floatalloc(nx);
    sax  = sf_complexalloc2(nx,nth);
    rax  = sf_complexalloc2(nx,nth);
  
    /* TAPER INIT */	
    wem2d_iso_rwetaper_init();
  
    /* FOURIER definitions */
    nk = nx;
    dk = 2.f*SF_PI/(nx*dx);
    ok = -nx/2.f*dk;
  
    /* Initialize Fourier transform */
    fft1_axis2_init(nx,nx);

    /* Initialize kernel */
    wem2d_iso_ker_init(nth,nx,dx,dz);
  
    /* Get minimum value for velocity filtering */
    float m=SF_HUGE;
    ZXLOOP( if (vel[iz][ix] < m) m = vel[iz][ix]; );
    ZLOOP( vmin[iz]=m; );
    
}

/*-----------------------------------------------------------------*/
/* WEMIG propgation */
void wem2d_iso_wemig(bool adj,
		     bool add,
		     float **vel, 
		     sf_complex **swf,
		     sf_complex **rwf,
		     float ***xig)
/*< Imaging kernel >*/
{
    int iw,ix,iz,ith,ih=0,rind,sind;
    int id=0;
    float sarg,rarg;
    sarg = 1.f;
    rarg =-1.f;
   
    if (adj) {

	/* Frequency Loop */  
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) private(iw,id,iz,ix,ih,rind,sind)
#endif
	for (iw = 0; iw < nw; iw++){
#ifdef _OPENMP
	    id = omp_get_thread_num();
#endif

	    fprintf(stderr,"\b\b\b\b\b\b\b\b\b\b\b\b%d of %d",iw+1,nw); 
    
	    /* Frequency for thread */
	    ww[id] = iw*dw+ow;
      
	    /* Input SWF and RWF */
	    XLOOP( sax[id][ix] = swf[iw][ix]; );
	    XLOOP( rax[id][ix] = rwf[iw][ix]; );
      
	    /* Depth Loop */
	    for (iz = 0; iz < nz; iz++) {
	
		/* WAVEFIELD TAPER */
		wem2d_iso_rwetaper(rax,id);
		wem2d_iso_rwetaper(sax,id);
	
		/* Step wavefields */
		wem2d_iso_shot_ker_onestep(id,iz,rarg,rax,vel,ww);
		wem2d_iso_shot_ker_onestep(id,iz,sarg,sax,vel,ww);
        
		/* High-angle filter FFT to kx */
		fft1_axis2((kiss_fft_cpx **) sax,id);
		fft1_axis2((kiss_fft_cpx **) rax,id);
	
		/* Fourier domain filtering */
		wem2d_iso_phs_correction(id,iz,rax,ww,vmin,rarg);
		wem2d_iso_phs_correction(id,iz,sax,ww,vmin,sarg);
	
		/* High-angle filter FFT to kx */
		ifft1_axis2((kiss_fft_cpx **) rax,id);
		ifft1_axis2((kiss_fft_cpx **) sax,id);
	
		/* WAVEFIELD TAPER */
		wem2d_iso_rwetaper(rax,id);
		wem2d_iso_rwetaper(sax,id);
	
		/* XIG Imaging Condition */
		for (ix=0; ix<nx; ix++) {
		    for (ih=0; ih<nh; ih++) {
			sind = SF_MAX(SF_MIN(ix+(nh-1)/2-ih,nx-1),0);
			rind = SF_MAX(SF_MIN(ix-(nh-1)/2+ih,nx-1),0);
			txig[iz][id][ix][ih] += crealf( rax[id][rind]*conjf(sax[id][sind]));
		    }
		}
	
	    } /* END DEPTH LOOP */
      
	    /* OUTPUT CURRENT SWF AND RWF */
	    XLOOP( swf[iw][ix] = sax[id][ix]; );
	    XLOOP( rwf[iw][ix] = rax[id][ix]; );
	
	} /* END FREQUENCIES */
    
	/* Sum image over all frequencies */
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) private(iz,ith,ix,ih)
#endif
	ZNXHLOOP( xig[iz][ix][ih] += txig[iz][ith][ix][ih]; );

    } else { /* MODELLING */
    
	/* Frequency Loop */  
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) private(iw,id,iz,ix,ih,rind,sind)
#endif
	for (iw = 0; iw < nw; iw++){
#ifdef _OPENMP
	    id = omp_get_thread_num();
#endif
      
	    fprintf(stderr,"\b\b\b\b\b\b\b\b\b\b\b\b%d of %d",iw+1,nw); 
      
	    /* Frequency for thread */
	    ww[id] = iw*dw+ow;
     
	    XLOOP( sax[id][ix] = swf[iw][ix]; );
	    XLOOP( rax[id][ix] = rwf[iw][ix]; );

	    for (iz = nz-1; iz > -1; iz--) {

		/* Spray xig out to all frequencies */
		for (ih = 0; ih < nh; ih++)
			XLOOP(txig[iz][id][ix][ih] = xig[iz][ix][ih];);

		/* Adjoint Imaging Condition */
		//for (ix=0; ix<nx; ix++) {
		//  for (ih=0; ih<nh; ih++) {
		//    sind = SF_MAX(SF_MIN(ix-(nh-1)/2+ih,nx-1),0);
		//    sax[id][sind] += txig[iz][id][ix][ih] ;// * (rax[id][sind]); 
		//  }
		//}

		for (ix=0; ix<nx; ix++) {
		    for (ih=0; ih<nh; ih++) {
			sind = SF_MAX(SF_MIN(ix-(nh-1)/2+ih,nx-1),0);
			rind = SF_MAX(SF_MIN(ix+(nh-1)/2-ih,nx-1),0);
			rax[id][rind] += txig[iz][id][ix][ih] * (sax[id][rind]) ;
		    }
		}

		/* WAVEFIELD TAPER */
		wem2d_iso_rwetaper(sax,id);
		wem2d_iso_rwetaper(rax,id); 

		/* High-angle filter FFT to kx */
		ifft1_axis2((kiss_fft_cpx **) sax,id);
		ifft1_axis2((kiss_fft_cpx **) rax,id);

		/* Fourier domain filtering */
		wem2d_iso_phs_correction(id,iz,sax,ww,vmin,rarg);
		wem2d_iso_phs_correction(id,iz,rax,ww,vmin,sarg);

		/* High-angle filter FFT to kx */
		fft1_axis2((kiss_fft_cpx **) rax,id);
		fft1_axis2((kiss_fft_cpx **) sax,id);

		/* Step wavefields */
		wem2d_iso_shot_ker_onestep(id,iz,rarg,sax,vel,ww);
		wem2d_iso_shot_ker_onestep(id,iz,sarg,rax,vel,ww);

		/* WAVEFIELD TAPER */
		wem2d_iso_rwetaper(rax,id);
		wem2d_iso_rwetaper(sax,id);
	    } /* END DEPTH LOOP */

	    /* Get output wavefields */
	    XLOOP( rwf[iw][ix] = rax[id][ix]; );
	    XLOOP( swf[iw][ix] = sax[id][ix]; );
    
	} /* End w LOOP */

    } /* End FOR/ADJ LOOP */

    /* Clean up memory */
    wem2d_iso_shot_ker_release();
  
}

/*-----------------------------------------------------------------*/

/* Set up taper */
static void wem2d_iso_rwetaper_init(void) 
{
    int j1,ix;
  
    /* Initialize taper */
    XLOOP( tap[ix] = 1.; );
	
    /* Left boundary */
    for (ix = 0; ix < nxtap; ix++){
	j1 = abs(nxtap-ix-1);
	tap[ix] = cos ( SF_PI/2.f*j1/(nxtap-1.f) );
    }
  
    /* Right boundary */
    for (ix = 0; ix < nxtap; ix++){
	j1 = abs(nxtap-ix-1);
	tap[nx-ix-1] = cos (SF_PI/2.f*j1/(nxtap-1.f) );
    }	

    tap[0]=0.f; tap[nx-1]=0.f;
  
}

/*-----------------------------------------------------------------*/

/* Taper wavefield */
static void wem2d_iso_rwetaper(sf_complex **dax, int id)
{
    int ixx;

    for (ixx=0; ixx < nx; ixx++) {
	dax[id][ixx] = dax[id][ixx] * tap[ixx];
    }
}

/*-----------------------------------------------------------------*/


void wem2d_iso_close() 
/*< Free memory >*/
{
    free (tap);
    free (vmin);
    free (*sax); free (sax);
    free (*rax); free (rax);
    free (***txig); free(*txig); free(txig);	
  
}
