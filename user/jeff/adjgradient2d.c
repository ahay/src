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
#include "adjgradient2d.h"
/*^*/

static void adjgradient2d_rwetaper_init(void);
static void adjgradient2d_rwetaper(sf_complex **dax, int id);

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
static float *vmin; /* Min velocity as function of depth */
static float *tap; /* Boundary taper */
static sf_complex **sax; /* Single w SWF (OMP) */
static sf_complex **rax; /* Single w RWF (OMP) */
static sf_complex **saxadj; /* Single w adjoint SWF (OMP) */
static sf_complex **raxadj; /* Single w adjoint RWF (OMP) */
static float ***tgrd; /* gradient (OMP) */
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
void adjgradient2d_init(int  nx_in,  int nz_in,   int nw_in, int nh_in,
			float ox_in, float oz_in, float ow_in, float oh_in,
			float dx_in, float dz_in, float dw_in, float dh_in,
			int nxtap_in,float **vel)
/*< initialize >*/
{
    int iz,ix, ith;

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
    vmin= sf_floatalloc(nz);
    ww  = sf_floatalloc(nth);
    tap = sf_floatalloc(nx);
    sax = sf_complexalloc2(nx,nth);
    rax = sf_complexalloc2(nx,nth);
    saxadj = sf_complexalloc2(nx,nth);
    raxadj = sf_complexalloc2(nx,nth);
    tgrd = sf_floatalloc3(nx,nth,nz);

    for (iz=0; iz<nz; iz++) {
	for (ith=0; ith<nth; ith++) {
	    for (ix=0; ix<nx; ix++) {
		tgrd[iz][ith][ix]=0.f;
	    }
	}
    }

    /* TAPER INIT */	
    adjgradient2d_rwetaper_init();
  
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
void adjgradient2d_wemig( float **vel, 
			  float ***xig,
			  sf_complex **swf,
			  sf_complex **rwf,
			  float **grd)
/*< Imaging kernel >*/
{
    int iw,ix,iz,ih,ith,pind,mind,id=0;
    float  caus= 1.f;
    float acaus=-1.f;

    /*  if (down) { /\* DOWNWARD CONTINUATION *\/ */
  
    /*     /\* Frequency Loop *\/   */
    /* #ifdef _OPENMP */
    /* #pragma omp parallel for schedule(dynamic) private(iw,id,iz,ix) */
    /* #endif */
    /*     for (iw = 0; iw < nw; iw++){ */
    /* #ifdef _OPENMP */
    /*       id = omp_get_thread_num(); */
    /* #endif */
  
    /*       fprintf(stderr,"\b\b\b\b\b\b\b\b\b\b\b\b%d of %d",iw+1,nw);  */
      
    /*       /\* Frequency for thread *\/ */
    /*      ww[id] = iw*dw+ow; */
        
    /*       /\* Inject wavefield at surface *\/ */
    /*       XLOOP( dax[id][ix] = wfi[iw][0][ix]; ); */
   
    /*       /\*  Onput wavefield at surface *\/ */
    /*   XLOOP( wfo[iw][0][ix] = dax[id][ix]; ); */
  
    /*   /\* Depth Loop *\/ */
    /*   for (iz = 1; iz < nz; iz++) { */
  
    /* 	/\* WAVEFIELD TAPER *\/ */
    /* 	adjwfld2d_iso_rwetaper(dax,id); */
  
    /* 	/\* Step wavefields *\/ */
    /* 	wem2d_iso_shot_ker_onestep(id,iz,weisign,dax,vel,ww); */
  
    /* 	/\* High-angle filter FFT to kx *\/ */
    /* 	fft1_axis2(dax,id); */
  
    /* 	/\* Fourier domain filtering *\/ */
    /* 	wem2d_iso_phs_correction(id,iz,dax,ww,vmin,weisign); */
  
    /* 	/\* High-angle filter FFT to kx *\/ */
    /* 	ifft1_axis2(dax,id); */
  
    /* 	/\* WAVEFIELD TAPER *\/ */
    /* 	adjwfld2d_iso_rwetaper(dax,id); */
  
    /* 	/\* Inject at depth *\/ */
    /* 	XLOOP( dax[id][ix] += wfi[iw][iz][ix]; ); */
  
    /* 	/\* Output wavefield *\/ */
    /* 	XLOOP( wfo[iw][iz][ix] = dax[id][ix]; ); */
  
    /*   } /\* END IZ *\/ */
    /* } /\* END IW *\/ */
  
    /* } else { /\* UPWARD CONTINUATION *\/ */
  
    /* Frequency Loop */  
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) private(iw,id,iz,ix,ih,pind,mind) \
    shared (sax,rax,saxadj,raxadj,ww)
#endif
    for (iw = 0; iw < nw; iw++){
#ifdef _OPENMP
	id = omp_get_thread_num();
#endif
      
	fprintf(stderr,"\b\b\b\b\b\b\b\b\b\b\b\b%d of %d",iw+1,nw); 
      
	/* Frequency for thread */
	ww[id] = iw*dw+ow;
      
	/* Inject WFs at max depth */
	XLOOP( rax[id][ix] = rwf[iw][ix]; );
	XLOOP( sax[id][ix] = swf[iw][ix]; );
     
	XLOOP( saxadj[id][ix] = 0.f; );
	XLOOP( raxadj[id][ix] = 0.f; );

	/* Generate adjoint SWF at nz-1 */
	for (ix=0; ix<nx; ix++) {
	    for (ih=0; ih<nh; ih++) {
		mind = SF_MAX( SF_MIN (ix-(nh-1)/2+ih,nx-1),0);
		pind = SF_MAX( SF_MIN (ix+(nh-1)/2-ih,nx-1),0);
		saxadj[id][mind] += xig[nz-1][ix][ih]*rax[id][pind];
	    }
	}

	/* Generate adjoint RWF at nz-1 */
	for (ix=0; ix<nx; ix++) {
	    for (ih=0; ih<nh; ih++) {
		mind = SF_MAX( SF_MIN (ix-(nh-1)/2+ih,nx-1),0);
		pind = SF_MAX( SF_MIN (ix+(nh-1)/2-ih,nx-1),0);
		raxadj[id][pind] += xig[nz-1][ix][ih]*sax[id][mind] ;
	    }
	}

	/* Compute gradient */
	XLOOP( tgrd[nz-1][id][ix] += \
	       -cimagf(ww[id]*ww[id]*(				\
			   saxadj[id][ix]*conjf(sax[id][ix])+	\
			   rax[id][ix]*conjf(raxadj[id][ix]) )); );
			       
	for (iz = nz-2; iz > -1; iz--) {

	    /* Step wavefields */
	    wem2d_iso_shot_ker_onestep(id,iz,acaus,sax   ,vel,ww);
	    wem2d_iso_shot_ker_onestep(id,iz, caus,rax   ,vel,ww);
	    wem2d_iso_shot_ker_onestep(id,iz,acaus,saxadj,vel,ww);
	    wem2d_iso_shot_ker_onestep(id,iz, caus,raxadj,vel,ww);

	    /* High-angle filter FFT to kx */
	    fft1_axis2((kiss_fft_cpx **) sax,   id);
	    fft1_axis2((kiss_fft_cpx **) rax,   id);
	    fft1_axis2((kiss_fft_cpx **) saxadj,id);
	    fft1_axis2((kiss_fft_cpx **) raxadj,id);

	    /* Fourier domain filtering */
	    wem2d_iso_phs_correction(id,iz,sax   ,ww,vmin,acaus);
	    wem2d_iso_phs_correction(id,iz,rax   ,ww,vmin, caus);
	    wem2d_iso_phs_correction(id,iz,saxadj,ww,vmin,acaus);
	    wem2d_iso_phs_correction(id,iz,raxadj,ww,vmin, caus);

	    /* High-angle filter FFT to kx */
	    ifft1_axis2((kiss_fft_cpx **) sax,   id);
	    ifft1_axis2((kiss_fft_cpx **) rax,   id);
	    ifft1_axis2((kiss_fft_cpx **) saxadj,id);
	    ifft1_axis2((kiss_fft_cpx **) raxadj,id);

	    /* WAVEFIELD TAPER */
	    adjgradient2d_rwetaper(sax,   id); 
	    adjgradient2d_rwetaper(rax,   id);
	    adjgradient2d_rwetaper(saxadj,id);
	    adjgradient2d_rwetaper(raxadj,id);

	    /* Generate adjoint SWF */
	    for (ix=0; ix<nx; ix++) {
		for (ih=0; ih<nh; ih++) {
		    mind = SF_MAX( SF_MIN (ix-(nh-1)/2+ih,nx-1), 0);
		    pind = SF_MAX( SF_MIN (ix+(nh-1)/2-ih,nx-1), 0);
		    saxadj[id][mind] += xig[iz][ix][ih]*rax[id][pind];
		}
	    }
	
	    /* Generate adjoint RWF */
	    for (ix=0; ix<nx; ix++) {
		for (ih=0; ih<nh; ih++) {
		    mind = SF_MAX( SF_MIN (ix-(nh-1)/2+ih,nx-1), 0);
		    pind = SF_MAX( SF_MIN (ix+(nh-1)/2-ih,nx-1), 0);
		    raxadj[id][pind] += xig[iz][ix][ih]*sax[id][mind];
		}
	    }
	
	    /* Compute gradient */
	       		
	    XLOOP( tgrd[iz][id][ix]+=			    \
		   -cimagf(ww[id]*ww[id]*(			    \
			       saxadj[id][ix]*conj(sax[id][ix])+   \
			       rax[id][ix]*conj(raxadj[id][ix]) )););

	} /* END DEPTH LOOP */
   
    } /* End w LOOP */

    /* Sum over frequencies to obtain gradient */
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) private (iz,ith,ix)
#endif
    for (iz=0; iz<nz; iz++) {
	for (ith=0; ith<nth; ith++) {
	    for (ix=0; ix<nx; ix++) {
		grd[iz][ix] += tgrd[iz][ith][ix];
	    }
	}
    }


    //  } /* End FOR/ADJ LOOP */
      
    fprintf(stderr,"\n"); 

    /* Clean up memory */
    wem2d_iso_shot_ker_release();
  
}

/*-----------------------------------------------------------------*/

/* Set up taper */
static void adjgradient2d_rwetaper_init(void) 
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
static void adjgradient2d_rwetaper(sf_complex **dax, int id)
{
    int ixx;

    for (ixx=0; ixx < nx; ixx++) {
	dax[id][ixx] = dax[id][ixx] * tap[ixx];
    }
}

/*-----------------------------------------------------------------*/

void adjgradient2d_close(void)
/*< Free memory >*/ 
{
    free (tap);
    free (vmin);
    free (ww); 
    free (*sax);    free (sax   );
    free (*rax);    free (rax   );
    free (*saxadj); free (saxadj);
    free (*raxadj); free (raxadj);
    free (**tgrd);  free ( *tgrd); free(tgrd);
  
}
