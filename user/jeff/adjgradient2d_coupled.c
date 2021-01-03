/* Module for adjoint-state gradient calculation for image-domain WET */
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
#include "adjgradient2d_coupled.h"
/*^*/

static int nx; /* num inline */
static int nw; /* num frequencies */
static int nz; /* num depths */
static int nh; /* num horizontal XIG */

static float dx; /* inline sampling */
static float dz; /* depth sampling */
static float dw; /* frequency sampling */
static float dh; /* horizontal XIG sampling */

static float ox; /* inline offset */
static float oz; /* depth offset */
static float ow; /* frequency offset */
static float oh; /* Starting XIG offset */

static int nxtap; /* taper length */
static float *vmin1; /* Min velocity as function of depth */
static float *vmin2; /* Min velocity as function of depth */
static float *tap; /* Boundary taper */
static sf_complex **us1x; /* Single w SWF (OMP) */
static sf_complex **ur1x; /* Single w RWF (OMP) */
static sf_complex **as1x; /* Single w adjoint SWF (OMP) */
static sf_complex **ar1x; /* Single w adjoint RWF (OMP) */
static sf_complex **us2x; /* Single w SWF (OMP) */
static sf_complex **ur2x; /* Single w RWF (OMP) */
static sf_complex **as2x; /* Single w adjoint SWF (OMP) */
static sf_complex **ar2x; /* Single w adjoint RWF (OMP) */

static float ***tgrd1; /* gradient (OMP) */
static float ***tgrd2; /* gradient (OMP) */
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

/*-------------------------------------------------------------*/
/* Initialization routine */
void adjgradient2d_coupled_init(int  nx_in,  int nz_in,   int nw_in, int nh_in,
				float ox_in, float oz_in, float ow_in, float oh_in,
				float dx_in, float dz_in, float dw_in, float dh_in,
				int nxtap_in,float **vel1, float **vel2)
/*< initialize >*/
{
    int iz,ix,ith;
    float m;

    nx=nx_in; ox=ox_in; dx=dx_in;
    nz=nz_in; oz=oz_in; dz=dz_in; 
    nw=nw_in; ow=2.f*SF_PI*ow_in; dw=2.f*SF_PI*dw_in;
    nh=nh_in; dh=dh_in; oh=oh_in;

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
    vmin1= sf_floatalloc(nz);
    vmin2= sf_floatalloc(nz);
    ww  = sf_floatalloc(nth);
    tap = sf_floatalloc(nx);
    us1x = sf_complexalloc2(nx,nth);
    ur1x = sf_complexalloc2(nx,nth);
    us2x = sf_complexalloc2(nx,nth);
    ur2x = sf_complexalloc2(nx,nth);
    as1x = sf_complexalloc2(nx,nth);
    ar1x = sf_complexalloc2(nx,nth);
    as2x = sf_complexalloc2(nx,nth);
    ar2x = sf_complexalloc2(nx,nth);
    tgrd1 = sf_floatalloc3(nx,nth,nz);
    tgrd2 = sf_floatalloc3(nx,nth,nz);

    for (iz=0; iz<nz; iz++) {
	for (ith=0; ith<nth; ith++) {
	    for (ix=0; ix<nx; ix++) {
		tgrd1[iz][ith][ix]=0.f;
	    }
	}
    }

    for (iz=0; iz<nz; iz++) {
	for (ith=0; ith<nth; ith++) {
	    for (ix=0; ix<nx; ix++) {
		tgrd2[iz][ith][ix]=0.f;
	    }
	}
    }

    /* TAPER INIT */	
    adjgradient2d_coupled_rwetaper_init();
  
    /* FOURIER definitions */
    nk = nx;
    dk = 2.f*SF_PI/(nx*dx);
    ok = -nx/2.f*dk;
  
    /* Initialize Fourier transform */
    fft1_axis2_init(nx,nx);

    /* Initialize kernel */
    wem2d_iso_ker_init(nth,nx,dx,dz);
  
    /* Get minimum value for velocity filtering */
    for (iz=0; iz<nz; iz++) {
	m = SF_HUGE;  
	for (ix=0; ix<nx; ix++) {
	    if (vel1[iz][ix] < m) m = vel1[iz][ix];
	}
	vmin1[iz]=m;
    }	
  
    for (iz=0; iz<nz; iz++) {
	m = SF_HUGE;  
	for (ix=0; ix<nx; ix++) {
	    if (vel2[iz][ix] < m) m = vel2[iz][ix];
	}
	vmin2[iz]=m;
    }	 
 
 
}

/*-----------------------------------------------------------------*/
/* WEMIG propgation */
void adjgradient2d_coupled_wemig( float **vel1, 
				  float ***xig1,
				  sf_complex **us1,
				  sf_complex **ur1,
				  float **grd1,
				  float **vel2, 
				  float ***xig2,
				  sf_complex **us2,
				  sf_complex **ur2,
				  float **grd2)
/*< Imaging kernel >*/
{
    int iw,ix,iz,ih,ith,pind,mind,id=0;
    float  caus= 1.f;
    float acaus=-1.f;

    /* Frequency Loop */  
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) private(iw,id,iz,ix,ih,pind,mind) \
    shared (us1x,ur1x,as1x,ar1x,us2x,ur2x,as2x,ar2x,ww,tgrd1,tgrd2)
#endif
    for (iw = 0; iw < nw; iw++){
#ifdef _OPENMP
	id = omp_get_thread_num();
#endif
      
	fprintf(stderr,"\b\b\b\b\b\b\b\b\b\b\b\b%d of %d",iw+1,nw); 
      
	/* Frequency for thread */
	ww[id] = iw*dw+ow;
      
	/* Inject WFs at max depth */
	XLOOP( ur1x[id][ix] = ur1[iw][ix]; );
	XLOOP( us1x[id][ix] = us1[iw][ix]; );
	XLOOP( ur2x[id][ix] = ur2[iw][ix]; );
	XLOOP( us2x[id][ix] = us2[iw][ix]; );
     
	XLOOP( as1x[id][ix] = 0.f; );
	XLOOP( ar1x[id][ix] = 0.f; );
	XLOOP( as2x[id][ix] = 0.f; );
	XLOOP( ar2x[id][ix] = 0.f; );

	/* Generate adjoint SWF at nz-1 */
	for (ix=0; ix<nx; ix++) {
	    for (ih=0; ih<nh; ih++) {
		mind = SF_MAX( SF_MIN (ix-(nh-1)/2+ih,nx-1),0);
		pind = SF_MAX( SF_MIN (ix+(nh-1)/2-ih,nx-1),0);
		as1x[id][mind] += xig1[nz-1][ix][ih]*ur1x[id][pind];
		as2x[id][mind] += xig2[nz-1][ix][ih]*ur2x[id][pind];
	    }
	}

	/* Generate adjoint RWF at nz-1 */
	for (ix=0; ix<nx; ix++) {
	    for (ih=0; ih<nh; ih++) {
		mind = SF_MAX( SF_MIN (ix-(nh-1)/2+ih,nx-1),0);
		pind = SF_MAX( SF_MIN (ix+(nh-1)/2-ih,nx-1),0);
		ar1x[id][pind] += xig1[nz-1][ix][ih]*us1x[id][mind] ;
		ar2x[id][pind] += xig2[nz-1][ix][ih]*us2x[id][mind] ;
	    }
	}

	/* Compute gradients */
/*
  XLOOP( tgrd1[nz-1][id][ix] += \
  -cimagf(ww[id]*ww[id]*(				\
  as1x[id][ix]*conjf(us1x[id][ix])+	\
  ur1x[id][ix]*conjf(ar1x[id][ix]) )); );
  XLOOP( tgrd2[nz-1][id][ix] += \
  -cimagf(ww[id]*ww[id]*(				\
  as2x[id][ix]*conjf(us2x[id][ix])+	\
  ur2x[id][ix]*conjf(ar2x[id][ix]) )); );
*/
	for (ix=0;ix<nx;ix++){
	    tgrd1[nz-1][id][ix]+= \
		-1.f/ww[id]*cimagf(as1x[id][ix]*conj(us1x[id][ix])+ur1x[id][ix]*conj(ar1x[id][ix]));}

	for (ix=0;ix<nx;ix++){ 
	    tgrd2[nz-1][id][ix]+= \
		-1.f/ww[id]*cimagf(as2x[id][ix]*conj(us2x[id][ix])+ur2x[id][ix]*conj(ar2x[id][ix]));}
	
			    
	for (iz = nz-2; iz > -1; iz--) {

	    /* Step wavefields */
	    wem2d_iso_shot_ker_onestep(id,iz,acaus,us1x,vel1,ww);
	    wem2d_iso_shot_ker_onestep(id,iz, caus,ur1x,vel1,ww);
	    wem2d_iso_shot_ker_onestep(id,iz,acaus,as1x,vel1,ww);
	    wem2d_iso_shot_ker_onestep(id,iz, caus,ar1x,vel1,ww);
	    wem2d_iso_shot_ker_onestep(id,iz,acaus,us2x,vel2,ww);
	    wem2d_iso_shot_ker_onestep(id,iz, caus,ur2x,vel2,ww);
	    wem2d_iso_shot_ker_onestep(id,iz,acaus,as2x,vel2,ww);
	    wem2d_iso_shot_ker_onestep(id,iz, caus,ar2x,vel2,ww);

	    /* High-angle filter FFT to kx */
	    fft1_axis2((kiss_fft_cpx **)us1x,id);
	    fft1_axis2((kiss_fft_cpx **)ur1x,id);
	    fft1_axis2((kiss_fft_cpx **)as1x,id);
	    fft1_axis2((kiss_fft_cpx **)ar1x,id);
	    fft1_axis2((kiss_fft_cpx **)us2x,id);
	    fft1_axis2((kiss_fft_cpx **)ur2x,id);
	    fft1_axis2((kiss_fft_cpx **)as2x,id);
	    fft1_axis2((kiss_fft_cpx **)ar2x,id);
		
	    /* Fourier domain filtering */
	    wem2d_iso_phs_correction(id,iz,us1x,ww,vmin1,acaus);
	    wem2d_iso_phs_correction(id,iz,ur1x,ww,vmin1, caus);
	    wem2d_iso_phs_correction(id,iz,as1x,ww,vmin1,acaus);
	    wem2d_iso_phs_correction(id,iz,ar1x,ww,vmin1, caus);
	    wem2d_iso_phs_correction(id,iz,us2x,ww,vmin2,acaus);
	    wem2d_iso_phs_correction(id,iz,ur2x,ww,vmin2, caus);
	    wem2d_iso_phs_correction(id,iz,as2x,ww,vmin2,acaus);
	    wem2d_iso_phs_correction(id,iz,ar2x,ww,vmin2, caus);

	    /* High-angle filter FFT to kx */
	    ifft1_axis2((kiss_fft_cpx **)us1x,id);
	    ifft1_axis2((kiss_fft_cpx **)ur1x,id);
	    ifft1_axis2((kiss_fft_cpx **)as1x,id);
	    ifft1_axis2((kiss_fft_cpx **)ar1x,id);
	    ifft1_axis2((kiss_fft_cpx **)us2x,id);
	    ifft1_axis2((kiss_fft_cpx **)ur2x,id);
	    ifft1_axis2((kiss_fft_cpx **)as2x,id);
	    ifft1_axis2((kiss_fft_cpx **)ar2x,id);
		
	    /* WAVEFIELD TAPER */
	    adjgradient2d_coupled_rwetaper(us1x,id); 
	    adjgradient2d_coupled_rwetaper(ur1x,id);
	    adjgradient2d_coupled_rwetaper(as1x,id);
	    adjgradient2d_coupled_rwetaper(ar1x,id);
	    adjgradient2d_coupled_rwetaper(us2x,id); 
	    adjgradient2d_coupled_rwetaper(ur2x,id);
	    adjgradient2d_coupled_rwetaper(as2x,id);
	    adjgradient2d_coupled_rwetaper(ar2x,id);

	    /* Generate adjoint SWF */
	    for (ix=0; ix<nx; ix++) {
		for (ih=0; ih<nh; ih++) {
		    mind = SF_MAX( SF_MIN (ix-(nh-1)/2+ih,nx-1), 0);
		    pind = SF_MAX( SF_MIN (ix+(nh-1)/2-ih,nx-1), 0);
		    as1x[id][mind] += xig1[iz][ix][ih]*ur1x[id][pind];
		    as2x[id][mind] += xig2[iz][ix][ih]*ur2x[id][pind];
		}
	    }	
	
	    /* Generate adjoint RWF */
	    for (ix=0; ix<nx; ix++) {
		for (ih=0; ih<nh; ih++) {
		    mind = SF_MAX( SF_MIN (ix-(nh-1)/2+ih,nx-1), 0);
		    pind = SF_MAX( SF_MIN (ix+(nh-1)/2-ih,nx-1), 0);
		    ar1x[id][pind] += xig1[iz][ix][ih]*us1x[id][mind];
		    ar2x[id][pind] += xig2[iz][ix][ih]*us2x[id][mind];
	    	}
	    }
	
	    for (ix=0;ix<nx;ix++){
		tgrd1[iz][id][ix]+= \
		    -1.f/ww[id]*cimagf(as1x[id][ix]*conj(us1x[id][ix])+ur1x[id][ix]*conj(ar1x[id][ix]));}
	    for (ix=0;ix<nx;ix++){ 
		tgrd2[iz][id][ix]+= \
		    -1.f/ww[id]*cimagf(as2x[id][ix]*conj(us2x[id][ix])+ur2x[id][ix]*conj(ar2x[id][ix]));}			
		
	    /* Compute gradient */
/*
  XLOOP( tgrd1[iz][id][ix]+=			    \
  -cimagf(ww[id]*ww[id]*(			    \
  as1x[id][ix]*conj(us1x[id][ix])+   \
  ur1x[id][ix]*conj(ar1x[id][ix]) )););
			
  XLOOP( tgrd2[iz][id][ix]+=			    \
  -cimagf(ww[id]*ww[id]*(			    \
  as2x[id][ix]*conj(us2x[id][ix])+   \
  ur2x[id][ix]*conj(ar2x[id][ix]) )););
*/
	} /* END DEPTH LOOP */
   
    } /* End w LOOP */

    /* Sum over frequencies to obtain gradient */
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) private (iz,ith,ix)
#endif
    for (iz=0; iz<nz; iz++) {
    	for (ith=0; ith<nth; ith++) {
	    for (ix=0; ix<nx; ix++) {
		grd1[iz][ix] += tgrd1[iz][ith][ix];
	    }
      	}
    }
    
    
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) private (iz,ith,ix)
#endif
    for (iz=0; iz<nz; iz++) {
    	for (ith=0; ith<nth; ith++) {
	    for (ix=0; ix<nx; ix++) {
		grd2[iz][ix] += tgrd2[iz][ith][ix];
	    }
      	}
    }
      
    fprintf(stderr,"\n"); 

    /* Clean up memory */
    wem2d_iso_shot_ker_release();
  
}

/*-----------------------------------------------------------------*/


void adjgradient2d_coupled_rwetaper_init() 
/*< Set up taper >*/
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

void adjgradient2d_coupled_rwetaper(sf_complex **dax, int id)
/*< Taper wavefield >*/
{
    int ixx;

    for (ixx=0; ixx < nx; ixx++) {
	dax[id][ixx] = dax[id][ixx] * tap[ixx];
    }
}

/*-----------------------------------------------------------------*/

void adjgradient2d_coupled_close() 
/*< Free memory >*/
{
    free (tap);
    free (vmin1);free(vmin2);
    free (ww); 
    free (*us1x); free (us1x);
    free (*ur1x); free (ur1x);
    free (*as1x); free (as1x);
    free (*ar1x); free (ar1x);
    free (**tgrd1); free(*tgrd1); free(tgrd1);
    free (*us2x); free (us2x);
    free (*ur2x); free (ur2x);
    free (*as2x); free (as2x);
    free (*ar2x); free (ar2x);
    free (**tgrd2); free(*tgrd2); free(tgrd2);  
}
