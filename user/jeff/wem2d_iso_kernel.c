/* KERNEL for 2D ISOTROPIC FD migration */
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
#include <math.h>
#include <rsf.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "wem2d_iso_kernel.h"
/*^*/

#define trick 0.11
#define a1 0.040315157
#define b1 0.873981642
#define a2 0.457289566 
#define b2 0.222691983

static int nx; /* num inline */
static float dx; /* inline sampling */
static float dz; /* depth sampling */
static int nth;

/* Fourier components */
static float fkx;
static float dkx;

/* FD Coeffs */
static float ax1;
static float bx1;
static float ax2;
static float bx2;

/* ARRAY DEFINITION*/
static sf_complex **dax; /* Single w wavefield (OMP) */
static sf_complex **rbx; /* Single w wavefield (OMP) */
static sf_complex **rcx; /* Single w wavefield (OMP) */
static sf_complex **ux;  /* Single w wavefield (OMP) */
static sf_complex **vx;  /* Single w wavefield (OMP) */
static float **cax; /* Single w wavefield (OMP) */
static float **xav; /* Single w wavefield (OMP) */
static float  *kkx; /* Single w wavefield (OMP) */

/*-------------------------------------------------------------*/
/* DEFINES */

#define XLOOP(a) for (ix=0; ix<nx; ix++) { {a} }
/*-----------------------------------------------------------------*/

void wem2d_iso_ker_init(	int nth_in,
				int nx_in,
				float dx_in,
				float dz_in)
/*< initialize >*/
{
    nx = nx_in;
    dx = dx_in;
    dz = dz_in;
    nth = nth_in;
    /*FOURIER*/
    fkx=-SF_PI/dx;
    dkx=2.f*SF_PI/(nx*dx);
	
    /*FD STUFF */
    ax1 = a1 * dz / (2.f*dx*dx);
    ax2 = a2 * dz / (2.f*dx*dx);
    bx1 = b1 / (dx*dx);
    bx2 = b2 / (dx*dx);
	
    dax = sf_complexalloc2(nx,nth);
    rbx = sf_complexalloc2(nx,nth);
    rcx = sf_complexalloc2(nx,nth);
    ux  = sf_complexalloc2(nx,nth);
    vx  = sf_complexalloc2(nx,nth);
    cax =   sf_floatalloc2(nx,nth);
    xav =   sf_floatalloc2(nx,nth);
    kkx =   sf_floatalloc     (nx);
	
    /* SET UP KX map */
    klinemap(kkx,nx,dx,0);
	
}

/*-----------------------------------------------------------------*/

/* SUBROUTINE TO STEP WAVEFIELD ONE SAMPLE IN DEPTH */
void wem2d_iso_shot_ker_onestep(int id, 
				int iz, 
				float signn,
				sf_complex **wld, 
				float **vel, 
				float *w)
/*< single stepper >*/
{
    int ix;
  
    /* FIND AVERAGE VELOCITY OVER FREQUENCY */
    XLOOP ( xav[id][ix]=vel[iz][ix]/w[id]; );
  
    /* SPLIT-STEP PROPAGATION */
    wem2d_iso_shot_ker_ssf(id,iz,wld,w,vel,signn);
  
    /* FIRST FFD STEP */
    wem2d_iso_shot_ker_fd(id,wld,xav,-signn*ax1,bx1);
  
    /* SECOND FFD STEP */
    wem2d_iso_shot_ker_fd(id,wld,xav,-signn*ax2,bx2);
  
}

/*-----------------------------------------------------------------*/

/* SUBROUTINE FOR SPLIT-STEP EXTRAPOLATION */
void wem2d_iso_shot_ker_ssf(int id, 
			    int iz, 
			    sf_complex **wld, 
			    float *w, 
			    float **vel, 
			    float signn)
/*< split step >*/
{
    int ix;
    float phsft;
    float tmp = -signn * dz * w[id];
    XLOOP( phsft = tmp / vel[iz][ix];
	   wld[id][ix]=wld[id][ix]*sf_cmplx( cos(phsft),sin(phsft) ); )
	}

/*-----------------------------------------------------------------*/

/* SUBROUTINE FOR FFD STEP */
void wem2d_iso_shot_ker_fd(int id, 
			   sf_complex **wld, 
			   float **ca, 
			   float ax, 
			   float bx)
/*< fourier FD step >*/
{
    int ix;

    XLOOP(  ux[id][ix] = wld[id][ix]; );
  
    XLOOP( cax[id][ix] =  ca[id][ix]; );
  
    /* TRIDIAGONAL SOLVER */
    wem2d_iso_fd(ux,vx,ax,bx,cax,dax,rbx,rcx,id);
  
    XLOOP( wld[id][ix] = vx[id][ix]; );

}

/*-----------------------------------------------------------------*/

/* FD CALLER ROUTINE */
void wem2d_iso_fd(sf_complex **u, 
		  sf_complex **v, 
		  float a, 
		  float b, 
		  float **ca, 
		  sf_complex **ra, 
		  sf_complex **rb, 
		  sf_complex **rc, 
		  int id)
/*< FD set up >*/
{
    int ix;
    sf_complex url, urr, ror, rol;
  
    /* CASE WHERE IX=0 */
    url   = 0.f;
    urr   = u[id][1];
    ror   = sf_cmplx( trick+b*ca[id][0]*ca[id][0] , a*ca[id][0]);
    rol   = conjf(ror);
    v [id][0]  = u[id][0]+ror*(url+urr-2.f*u[id][0]);
    ra[id][0] = 1.f - 2.f*rol;
    rb[id][0] = rol;
    rc[id][0] = 0.;
  
    /* NON-END CASE */
    for (ix=1; ix < nx-1; ix++){
	url      = u[id][ix-1];
	urr      = u[id][ix+1];
	ror      = sf_cmplx( trick+b*ca[id][ix]*ca[id][ix], a*ca[id][ix] );
	rol      = conjf(ror);
	v [id][ix]    = u[id][ix]+ror*(url+urr-2.f*u[id][ix]);
	ra[id][ix]   = 1.f -2.f*rol;
	rb[id][ix]   = rol;
	rc[id][ix-1] = rol;
    }
  
    /* CASE WHERE IX=NX-1 */
    url     = u[id][nx-2];
    urr     = 0.f;
    ror     = sf_cmplx( trick + b*ca[id][nx-1]*ca[id][nx-1] , a*ca[id][nx-1] );
    rol     = conjf(ror);
    v [id][nx-1] = u[id][nx-1] + ror*(url+urr-2.f*u[id][nx-1]);
    ra[id][nx-1] = 1.f -2.f*rol;
    rb[id][nx-1] = 0.f; 
    rc[id][nx-2] = rol;
  
    /* CALL TRIDIAGONAL SOLVER */
    wem2d_iso_thryj(ra,rb,rc,v,id);
  
}

/*-----------------------------------------------------------------*/

void wem2d_iso_thryj(sf_complex **a, 
		     sf_complex **b, 
		     sf_complex **c, 
		     sf_complex **v, 
		     int id)
/*< Tridiag solver >*/
{
    int ix;

    /* for i=0 */
    b[id][0] = b[id][0] / a[id][0]; 
    v[id][0] = v[id][0] / a[id][0]; 
  
    /* for i=1,nx-2 */
    for (ix=1; ix < nx-1; ix++){
	a[id][ix] = a[id][ix] - b[id][ix-1]*c[id][ix-1];
	v[id][ix] =(v[id][ix] - c[id][ix-1]*v[id][ix-1])/a[id][ix];
	b[id][ix] = b[id][ix] / a[id][ix];
    }
  
    /* for i=n-1 */
    a[id][nx-1] = a[id][nx-1] - b[id][nx-2]*c[id][nx-2];
    v[id][nx-1] =(v[id][nx-1] - c[id][nx-2]*v[id][nx-2])/a[id][nx-1];
  
    /* Back-substitution */
    for (ix=nx-2; ix > -1; ix--){
	v[id][ix] = v[id][ix] - b[id][ix]*v[id][ix+1];
    }   
  
}

/*-----------------------------------------------------------------*/

/* FOURIER-DOMAIN PHASE SHIFT */
void wem2d_iso_phs_correction(int id, 
			      int iz, 
			      sf_complex **dax, 
			      float *w, 
			      float *vmin, 
			      int darg)
/*< Phase shift correction >*/
{
    int ikx;

    float sr2, phsft, tx;
    float w_v0 = w[id] / vmin[iz];
  
    float wv2 = w_v0*w_v0;
    float avX = 1.f/w_v0;
  
    for (ikx=0; ikx < nx; ikx++){
	sr2 = kkx[ikx] * kkx[ikx] / wv2;	
    
	if (sr2 > 1.f ) {
	    dax[id][ikx]=0.f;
	} else {
	    tx = kkx[ikx]*kkx[ikx]*avX*avX;
	    phsft=-dz*w_v0*( sqrt(1.f-sr2)-(1.f-a1*tx/(1.f-b1*tx)-a2*tx/(1.f-b2*tx)));
	    dax[id][ikx] = dax[id][ikx]*sf_cmplx( cos(phsft) , darg*sin(phsft) );
	}
    }
}

/*-----------------------------------------------------------------*/

/* FOURIER FREQUENCY MAP */
/* NOT CENTERED */
void klinemap(float *kxmap, 
	      int n, 
	      float d, 
	      int shift)
/*< Generate wavenumber line >*/
{
    int i, m=n/2;
    float kscale=2.f*SF_PI/(n*d);
  
    for (i=0; i<m; i++){
	kxmap[i+shift]=kscale*(i-1.f);
    }
  
    for (i=m+1; i<n; i++){
	kxmap[i+shift]=kscale*(i-n-1.f);
    }
  
}

/*-----------------------------------------------------------------*/

void wem2d_iso_shot_ker_release()
/*< Free memory >*/
{
    free(*dax); free(dax);
    free(*rbx); free(rbx);
    free(*rcx); free(rcx);
    free(*ux);  free(ux);
    free(*vx);  free(vx);
    free(*cax); free(cax);
    free(*xav); free(xav);
    free(kkx);
}

