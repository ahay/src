/* 3-D linearized SSR */

/*
  Copyright (C) 2007 Colorado School of Mines
  
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
/*^*/

#include "ompfft.h"
#include "taper3.h"

#include "weutil.h"
/*^*/

/*------------------------------------------------------------*/

#define LOOP(a) for(iy=0;iy<cub->amy.n;iy++){	\
	        for(ix=0;ix<cub->amx.n;ix++){ \
		    {a} }}

#define KOOP(a) for(iy=0;iy<lsr->byy.n;iy++){ \
	        for(ix=0;ix<lsr->bxx.n;ix++){ \
		    {a} }}

#define KMAP(i,n) (i<n/2.) ? (i+n/2.) : (i-n/2.);
#define X2K(a,b,p) b.n=a.n+p; b.d=2.0*SF_PI/(b.n*a.d); b.o=(1==b.n)?0:-SF_PI/a.d;

/*static sf_axa az,axx,ayy;*/
/*static sf_axa    bxx,byy;*/

/*static float         **kk; */
/* wavenumber  */
/*static float         **kw; */
/* wavenumber weight */

/*static sf_complex **wk;*/
/*static sf_complex **wt;*/

static int   nsc;
static float csc[5];

/*------------------------------------------------------------*/
lsr3d lsr3_init( cub3d cub,
		 int px,
		 int py
    )
/*< initialize >*/
{
    int   ix, iy;
    float kx, ky;
    int   jx, jy;

    /*------------------------------------------------------------*/
    lsr3d lsr;
    lsr = (lsr3d) sf_alloc(1,sizeof(*lsr));

    X2K(cub->amx,lsr->bxx,px);
    X2K(cub->amy,lsr->byy,py);

    /* allocate K-domain storage */
    lsr->kk = sf_floatalloc2(lsr->bxx.n,lsr->byy.n);
    lsr->kw = sf_floatalloc2(lsr->bxx.n,lsr->byy.n);
    lsr->wk = sf_complexalloc3(lsr->bxx.n,lsr->byy.n,cub->ompnth);
    lsr->wt = sf_complexalloc3(lsr->bxx.n,lsr->byy.n,cub->ompnth);

    for (iy=0; iy<lsr->byy.n; iy++) {
	jy = KMAP(iy,lsr->byy.n);
	ky = lsr->byy.o + jy*lsr->byy.d;
	for (ix=0; ix<lsr->bxx.n; ix++) {
	    jx = KMAP(ix,lsr->bxx.n);
	    kx = lsr->bxx.o + jx*lsr->bxx.d;  
	    lsr->kk[iy][ix] = kx*kx+ky*ky;

	    lsr->kw[iy][ix] = 1.;
	}
    }    

    /* initialize FFT */
    lsr->f2d = ompfft2_init(cub,lsr->bxx.n,lsr->byy.n);

    /* square-root expansion coefficients */
    if(!sf_getint("nsc",&nsc)) nsc = 0;
    if(nsc>5) nsc=5;

    csc[0]= 1.;
    csc[1]= 1./  2.;
    csc[2]= 3./  8.;
    csc[3]= 5./ 16.;
    csc[4]=35./128.;

    return lsr;
}

/*------------------------------------------------------------*/

void lsr3_close(lsr3d lsr)
/*< free allocated storage >*/
{
    ompfft2_close(lsr->f2d);

    free(*lsr->kk); free(lsr->kk);
    free(*lsr->kw); free(lsr->kw);

    free(**lsr->wk); free(*lsr->wk); free(lsr->wk);
    free(**lsr->wt); free(*lsr->wt); free(lsr->wt);
}

/*------------------------------------------------------------*/

void kweight( cub3d  cub,
	      lsr3d  lsr,
	      float **bs, /* slowness */
	      float   wo  /* frequency */
)
/*< k-domain weight >*/
{
    int ix,iy,ii,nn;
    float *ss, smin, ko;
    
    /* find s min */
    nn = cub->amx.n*cub->amy.n;
    ss = sf_floatalloc(nn);
    ii=0;
    LOOP( ss[ii] = bs[iy][ix];
	  ii++; );
    smin = sf_quantile(0,nn,ss);
    free(ss);

    ko  = abs(wo) * smin;
    ko *= ko;

    KOOP(
	if( lsr->kk[iy][ix] < ko ) {
	    lsr->kw[iy][ix] = 1.;
	} else {
	    lsr->kw[iy][ix] = 0.;
	} );
}

/*------------------------------------------------------------*/
void lsr3_s2w(cub3d  cub,
	      lsr3d  lsr,
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
    int ompith=0;

    wo = cimagf(w);     /* real frequency */
    iwdz = sf_cmplx(0.,-2 * wo * cub->amz.d);

    /* 0th order term */
#ifdef SF_HAS_COMPLEX_H
    LOOP( pw[iy][ix] =
	  ps[iy][ix] * iwdz * lsr->bw[ompith][iy][ix]; );
#else
    LOOP( pw[iy][ix] =
	  sf_cmul(ps[iy][ix], sf_cmul(iwdz, lsr->bw[ompith][iy][ix])); );
#endif

    /* higher order terms */
    if(nsc>0) {
	kweight(cub,lsr,bs,wo);     /* k-domain weight */

	LOOP( lsr->wt[iy][ix][ompith] = pw[iy][ix]; );

	for( isc=1; isc<=nsc; isc++) {
	    KOOP( lsr->wk[iy][ix][ompith] = sf_cmplx(0.,0.); );
	    LOOP( lsr->wk[iy][ix][ompith] = lsr->wt[iy][ix][ompith]; );

	    ompfft2(false,(kiss_fft_cpx**) lsr->wk[ompith],ompith,lsr->f2d);

#ifdef SF_HAS_COMPLEX_H
	    KOOP( lsr->wk[iy][ix][ompith] *= 
		  lsr->kw[iy][ix] * pow(lsr->kk[iy][ix],isc); );
#else
	    KOOP( lsr->wk[iy][ix] = sf_crmul(lsr->wk[iy][ix], 
					lsr->kw[iy][ix] * pow(lsr->kk[iy][ix],isc)); );
#endif

	    ompfft2(true,(kiss_fft_cpx**) lsr->wk[ompith],ompith,lsr->f2d);

#ifdef SF_HAS_COMPLEX_H
	    LOOP( lsr->wk[iy][ix][ompith] *= pow(wo*bs[iy][ix],-2*isc); );
	    LOOP( pw[iy][ix] +=
		  lsr->wk[iy][ix][ompith] * csc[isc]; );
#else
	     LOOP( lsr->wk[iy][ix][ompith] = sf_crmul(lsr->wk[iy][ix][ompith],
					 pow(wo*bs[iy][ix],-2*isc)); );
	     LOOP( pw[iy][ix] = sf_cadd(pw[iy][ix],
					sf_crmul(lsr->wk[iy][ix],csc[isc])); );
#endif
	}
    }
}

/*------------------------------------------------------------*/
void lsr3_w2s(cub3d  cub,
	      lsr3d  lsr,
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
    int ompith=0;

    wo = cimagf(w);     /* real frequency */
    iwdz = sf_cmplx(0.,-2 * wo * cub->amz.d);

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
	kweight(cub,lsr,bs,wo);     /* k-domain weight */

	for( isc=1; isc<=nsc; isc++) {

	    KOOP( lsr->wk[iy][ix][ompith] = sf_cmplx(0.,0.); );
	    LOOP( lsr->wk[iy][ix][ompith] = pw[iy][ix]; );
#ifdef SF_HAS_COMPLEX_H
	    LOOP( lsr->wk[iy][ix][ompith]*= pow(wo*bs[iy][ix],-2*isc); );
#else
	    LOOP( lsr->wk[iy][ix][ompith] = sf_crmul(lsr->wk[iy][ix],
					pow(wo*bs[iy][ix],-2*isc)); );
#endif

	    ompfft2(true,(kiss_fft_cpx**) lsr->wk[ompith],ompith,lsr->f2d);

#ifdef SF_HAS_COMPLEX_H
	    KOOP( lsr->wk[iy][ix][ompith] *=
		  lsr->kw[iy][ix] * pow(lsr->kk[iy][ix],isc); );
#else
	    KOOP( lsr->wk[iy][ix][ompith] = sf_crmul(lsr->wk[iy][ix],
					lsr->kw[iy][ix] * pow(lsr->kk[iy][ix],isc)); );
#endif

	    ompfft2(false,(kiss_fft_cpx**) lsr->wk[ompith],ompith,lsr->f2d);

#ifdef SF_HAS_COMPLEX_H	    
	    LOOP( lsr->wk[iy][ix][ompith] *= iwdz * conjf( bw[iy][ix] ); );
	    
	    LOOP( ps[iy][ix] +=
		  lsr->wk[iy][ix][ompith] * csc[isc]; );
#else
	    LOOP( lsr->wk[iy][ix][ompith] = sf_cmul(lsr->wk[iy][ix],
				       sf_cmul(iwdz,conjf( bw[iy][ix] ))); );
	    
	    LOOP( ps[iy][ix] = sf_cadd(ps[iy][ix],
				       sf_crmul(lsr->wk[iy][ix],csc[isc])); );
#endif
	}
    }
}














