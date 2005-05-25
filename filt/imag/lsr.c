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

#include "fft2.h"

#define LOOP(a) for(iy=0;iy<ayy.n;iy++){ for(ix=0;ix<axx.n;ix++){ {a} }}
#define KOOP(a) for(iy=0;iy<byy.n;iy++){ for(ix=0;ix<bxx.n;ix++){ {a} }}

#define  KMAP(i,n) (i<n/2.) ? (i+n/2.) : (i-n/2.);
#define X2K(a,b,p) b.n=a.n+p; b.d=2.0*SF_PI/(b.n*a.d); b.o=(1==b.n)?0:-SF_PI/a.d;

static axa az,axx,ayy;
static axa    bxx,byy;

static float         **kk; /* wavenumber  */
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
    for (iy=0; iy<byy.n; iy++) {
	jy = KMAP(iy,byy.n);
	ky = byy.o + jy*byy.d;
	for (ix=0; ix<bxx.n; ix++) {
	    jx = KMAP(ix,bxx.n);
	    kx = bxx.o + jx*bxx.d;  
	    kk[iy][ix] = kx*kx+ky*ky;
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

}

/*------------------------------------------------------------*/


void lsr_one(
    float complex    w /* frequency */,
    complex float **bw /* background wavefield */,
    complex float **ps /* perturbation slowness */,
    complex float **pw /* perturbation wavefield */,
    float          *sm /* ref slo squared */
    )
/*< linear scattering operator >*/
{
    float complex w2,cc;
    int ix,iy;

    w2 = w*w;

    LOOP( wk[iy][ix] = bw[iy][ix] * ps[iy][ix]; );

    fft2(false,wk);

    KOOP(
	cc = w / csqrtf(1+ kk[iy][ix] / (w2*sm[0]) );
	wk[iy][ix] *= cc;
	);

    fft2(true,wk);

    LOOP( pw[iy][ix] = wk[iy][ix]; );

}

void lsr_tst(
    float complex    w /* frequency */,
    complex float **bw /* background wavefield */,
    complex float **ps /* perturbation slowness */,
    complex float **pw /* perturbation wavefield */,
    float          *sm /* ref slo squared */
    )
/*< linear scattering operator >*/
{
    int ix,iy;

    float wo;
/*    float ko;*/
/*    float co;*/
/*    float xx;*/

    wo = cimagf(w);
/*    ko = wo*wo * sm[0];*/

    LOOP( wk[iy][ix] =
	  bw[iy][ix] * 
	  ps[iy][ix] * I * wo * az.d; );

/*    fft2(false,wk);*/
/*    */
/*    KOOP( xx = kk[iy][ix]/ko;*/
/*	  if(xx<0.5) {*/
/*	      co = 1/sqrtf(1-xx);*/
/*	      co = 1 + 0.5 * xx;*/
/*	  } else {*/
/*	      co = 0;*/
/*	  }*/
/*	  wk[iy][ix] *= co;*/
/*	  fprintf(stderr,"%g %g %g %g\n",kk[iy][ix],ko,xx,co);*/
/*	);*/
/**/
/*    fft2(true,wk);*/

    LOOP( pw[iy][ix] = wk[iy][ix]; );
}

void lsr_exp(
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

    wo = cimagf(w);

    KOOP( wk[iy][ix]=0.;
	  wt[iy][ix]=0.;
	);

    LOOP( pw[iy][ix] = 0.;
	  wk[iy][ix] =
	  bw[iy][ix] * 
	  ps[iy][ix] * I * wo * az.d; 
	);

    fft2(false,wk);
    
/*    fprintf(stderr,"nsc=%d\n",nsc);*/
    for( isc=0; isc<nsc; isc++) {
	KOOP( wt[iy][ix] = 
	      wk[iy][ix] * pow(kk[iy][ix],isc); 
	    );
    
	fft2(true,wt);
    
	LOOP( pw[iy][ix] +=
	      wt[iy][ix] * 
	      csc[isc] / pow(wo*bs[iy][ix],2*isc);
	);
    }
}
