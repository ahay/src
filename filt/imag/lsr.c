/* 3-D linearized SSR */

/*
  Copyright (C) 2006 Colorado School of Mines
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

#define KMAP(i,n) (i<n/2.) ? (i+n/2.) : (i-n/2.);
#define X2K(a,b,p) b.n=a.n+p; b.d=2.0*SF_PI/(b.n*a.d); b.o=(1==b.n)?0:-SF_PI/a.d;

static sf_axa az,axx,ayy;
static sf_axa    bxx,byy;

static float         **kk; /* wavenumber  */
static float         **kw; /* wavenumber weight */

static sf_complex **wk;
static sf_complex **wt;

static int   nsc;
static float csc[5];

/*------------------------------------------------------------*/

void lsr_init( sf_axis az_,
	       sf_axis ax_,
	       sf_axis ay_,
	       int px,
	       int py
    )
/*< initialize >*/
{
    int   ix, iy;
    float kx, ky;
    int   jx, jy;

    az  = sf_nod(az_);
    axx = sf_nod(ax_);
    ayy = sf_nod(ay_);

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
    if(!sf_getint("nsc",&nsc)) nsc = 0;
    if(nsc>5) nsc=5;

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

    /* find s min */
    nn = axx.n*ayy.n;
    ss = sf_floatalloc(nn);
    ii=0;
    LOOP( ss[ii] = bs[iy][ix];
	  ii++; );
    smin = sf_quantile(0,nn,ss);
    free(ss);

    ko  = abs(wo) * smin;
    ko *= ko;

    KOOP(
	if( kk[iy][ix] < ko ) {
	    kw[iy][ix] = 1.;
	} else {
	    kw[iy][ix] = 0.;
	} );
}

/*------------------------------------------------------------*/
void lsr_s2w(
    sf_complex    w /* frequency */,
    sf_complex **bw /* background   wavefield */,
    float      **bs /* background   slowness  */,
    sf_complex **pw /* perturbation wavefield */,
    sf_complex **ps /* perturbation slowness  */
    )
/*< linear scattering operator (forward) >*/
{
    int ix,iy,isc;
    float wo;
    sf_complex iwdz;

    wo = cimagf(w);     /* real frequency */
    iwdz = sf_cmplx(0.,-2 * wo * az.d);

    /* 0th order term */
#ifdef SF_HAS_COMPLEX_H
    LOOP( pw[iy][ix] =
	  ps[iy][ix] * iwdz * bw[iy][ix]; );
#else
    LOOP( pw[iy][ix] =
	  sf_cmul(ps[iy][ix], sf_cmul(iwdz, bw[iy][ix])); );
#endif

    /* higher order terms */
    if(nsc>0) {
	kweight(bs,wo);     /* k-domain weight */

	LOOP( wt[iy][ix] = pw[iy][ix]; );

	for( isc=1; isc<=nsc; isc++) {
	    KOOP( wk[iy][ix] = sf_cmplx(0.,0.); );
	    LOOP( wk[iy][ix] = wt[iy][ix]; );

	    fft2(false,(kiss_fft_cpx**) wk);

#ifdef SF_HAS_COMPLEX_H
	    KOOP( wk[iy][ix] *= 
		  kw[iy][ix] * pow(kk[iy][ix],isc); );
#else
	    KOOP( wk[iy][ix] = sf_crmul(wk[iy][ix], 
					kw[iy][ix] * pow(kk[iy][ix],isc)); );
#endif

	    fft2(true,(kiss_fft_cpx**) wk);

#ifdef SF_HAS_COMPLEX_H
	    LOOP( wk[iy][ix] *= pow(wo*bs[iy][ix],-2*isc); );
	    LOOP( pw[iy][ix] +=
		  wk[iy][ix] * csc[isc]; );
#else
	     LOOP( wk[iy][ix] = sf_crmul(wk[iy][ix],
					 pow(wo*bs[iy][ix],-2*isc)); );
	     LOOP( pw[iy][ix] = sf_cadd(pw[iy][ix],
					sf_crmul(wk[iy][ix],csc[isc])); );
#endif
	}
    }
}

/*------------------------------------------------------------*/
void lsr_w2s(
    sf_complex    w /* frequency */,
    sf_complex **bw /* background   wavefield */,
    float      **bs /* background   slowness  */,
    sf_complex **pw /* perturbation wavefield */,
    sf_complex **ps /* perturbation slowness  */
    )
/*< linear scattering operator (adjoint) >*/
{
    int ix,iy,isc;
    float wo;
    sf_complex iwdz;

    wo = cimagf(w);     /* real frequency */
    iwdz = sf_cmplx(0.,-2 * wo * az.d);

    /* 0th order term */
#ifdef SF_HAS_COMPLEX_H
    LOOP( ps[iy][ix] =
	  pw[iy][ix] * iwdz * conjf( bw[iy][ix] ); );
#else
    LOOP( ps[iy][ix] =
	  sf_cmul(pw[iy][ix],sf_cmul(iwdz,conjf( bw[iy][ix] ))); );
#endif

    /* higher order terms */
    if(nsc>0) {
	kweight(bs,wo);     /* k-domain weight */

	for( isc=1; isc<=nsc; isc++) {
	    KOOP( wk[iy][ix] = sf_cmplx(0.,0.); );
	    LOOP( wk[iy][ix] = pw[iy][ix]; );
#ifdef SF_HAS_COMPLEX_H
	    LOOP( wk[iy][ix]*= pow(wo*bs[iy][ix],-2*isc); );
#else
	    LOOP( wk[iy][ix] = sf_crmul(wk[iy][ix],
					pow(wo*bs[iy][ix],-2*isc)); );
#endif

	    fft2(true,(kiss_fft_cpx**) wk);

#ifdef SF_HAS_COMPLEX_H
	    KOOP( wk[iy][ix] *=
		  kw[iy][ix] * pow(kk[iy][ix],isc); );
#else
	    KOOP( wk[iy][ix] = sf_crmul(wk[iy][ix],
					kw[iy][ix] * pow(kk[iy][ix],isc)); );
#endif

	    fft2(false,(kiss_fft_cpx**) wk);

#ifdef SF_HAS_COMPLEX_H	    
	    LOOP( wk[iy][ix] *= iwdz * conjf( bw[iy][ix] ); );
	    
	    LOOP( ps[iy][ix] +=
		  wk[iy][ix] * csc[isc]; );
#else
	    LOOP( wk[iy][ix] = sf_cmul(wk[iy][ix],
				       sf_cmul(iwdz,conjf( bw[iy][ix] ))); );
	    
	    LOOP( ps[iy][ix] = sf_cadd(ps[iy][ix],
				       sf_crmul(wk[iy][ix],csc[isc])); );
#endif
	}
    }
}














