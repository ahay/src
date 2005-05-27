/* 3-D linearized SSR */
/*
  Copyright (C) 2004 University of Texas at Austin
  
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

#include "fft2.h"

#define LOOP(a) for(iy=0;iy<ayy.n;iy++){ for(ix=0;ix<axx.n;ix++){ {a} }}
#define KOOP(a) for(iy=0;iy<byy.n;iy++){ for(ix=0;ix<bxx.n;ix++){ {a} }}

#define  KMAP(i,n) (i<n/2.) ? (i+n/2.) : (i-n/2.);
#define X2K(a,b,p) b.n=a.n+p; b.d=2.0*SF_PI/(b.n*a.d); b.o=(1==b.n)?0:-SF_PI/a.d;

static axa az,axx,ayy;
static axa    bxx,byy;

static float         **kk; /* wavenumber  */
static float         **kw; /* wavenumber weight */

static float complex **wk;
static float complex **wt;

static int   nsc;
static float csc[5];

/*------------------------------------------------------------*/

void lsr_init( axa az_,
	       axa ax_,
	       axa ay_,
	       int px,
	       int py
    )
/*< initialize >*/
{
    int   ix, iy;
    float kx, ky;
    int   jx, jy;

    az  = az_;
    axx = ax_;
    ayy = ay_;

    /* construct K-domain axes */
    X2K(axx,bxx,px);
    X2K(ayy,byy,py);
    fft2_init(bxx.n,byy.n);

    /* precompute wavenumbers */
    kk = sf_floatalloc2(bxx.n,byy.n);
    kw = sf_floatalloc2(bxx.n,byy.n);
    for (iy=0; iy<byy.n; iy++) {
	jy = KMAP(iy,byy.n);
	ky = byy.o + jy*byy.d;
	for (ix=0; ix<bxx.n; ix++) {
	    jx = KMAP(ix,bxx.n);
	    kx = bxx.o + jx*bxx.d;  
	    kk[iy][ix] = kx*kx+ky*ky;

	    kw[iy][ix] = 1.;
	}
    }    

    /* allocate K-domain storage */
    wk = sf_complexalloc2 (bxx.n,byy.n);
    wt = sf_complexalloc2 (bxx.n,byy.n);

    /* square-root expansion coefficients */
    nsc = 2;
    csc[0]= 1.;
    csc[1]= 1./  2.;
    csc[2]= 3./  8.;
    csc[3]= 5./ 16.;
    csc[4]=35./128.;
}

/*------------------------------------------------------------*/

void lsr_close(void)
/*< free allocated storage >*/
{
    free(*kk); free(kk);
    free(*kw); free(kw);

    free(*wk); free(wk);
    free(*wt); free(wt);
}

/*------------------------------------------------------------*/

void kweight( float **bs, /* slowness */
	      float   wo  /* frequency */
)
/*< k-domain weight >*/
{
    int ix,iy,ii,nn;
    float *ss, smin, ko;

    nn = axx.n*ayy.n;
    ss = sf_floatalloc(nn);
    ii=0;
    LOOP( ss[ii] = bs[iy][ix];
	  ii++; );
    smin = sf_quantile(0,nn,ss);
    free(ss);

    ko = wo * smin;

    KOOP( 
	if(kk[iy][ix] < ko) {
	    kw[iy][ix]=1.;
	} else {
	    kw[iy][ix]=0.;
	} );

}

/*------------------------------------------------------------*/

void lsr_s2w(
    complex float    w /* frequency */,
    complex float **bw /* background wavefield   */,
    float         **bs /* background slowness    */,
    complex float **pw /* perturbation wavefield */,
    complex float **ps /* perturbation slowness  */
    )
/*< linear scattering operator >*/
{
    int ix,iy;
    float wo;
    int  isc;
    float complex iwdz;

    wo = cimagf(w);     /* real frequency */
    iwdz = 2 * I * wo * az.d;

    kweight(bs,wo);     /* k-domain weight */

    KOOP( wk[iy][ix]=0.;
	  wt[iy][ix]=0.; );

    /* 0th order term */
    LOOP( pw[iy][ix] =
	  ( iwdz ) * 
	  bw[iy][ix] * 
	  ps[iy][ix]; );

    /* higher order terms */
    LOOP( wk[iy][ix] = pw[iy][ix]; );
    if(nsc>0) {
	fft2(false,wk);	
	for( isc=1; isc<=nsc; isc++) {
	    KOOP( wt[iy][ix] = 
		  kw[iy][ix] * 
		  wk[iy][ix] * pow(kk[iy][ix],isc); );

	    fft2(true,wt);
	    
	    LOOP( wt[iy][ix] /=  pow(wo*bs[iy][ix],2*isc); );
	    
	    LOOP( pw[iy][ix] +=
		  wt[iy][ix] * csc[isc]; );
	}
    } /* series summation */
}

void lsr_w2s(
    complex float    w /* frequency */,
    complex float **bw /* background   wavefield   */,
    float         **bs /* background   slowness    */,
    complex float **pw /* perturbation wavefield */,
    complex float **ps /* perturbation slowness  */
    )
/*< linear scattering operator >*/
{
    int ix,iy;
    float wo;
    float complex iwdz;

    wo = cimagf(w);     /* real frequency */
    iwdz = 2 * I * wo * az.d;

    kweight(bs,wo);     /* k-domain weight */

    /* 0th order term */
    LOOP( ps[iy][ix] =
	  conjf( bw[iy][ix] ) * 
	  conjf( iwdz ) * 
	  pw[iy][ix]; );

    /* higher order terms */
    
}
