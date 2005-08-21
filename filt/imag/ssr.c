/* 
 * 3-D SSR
 * pcs 2005
 */

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
/*^*/

#include "fft2.h"
#include "taper.h"

/*#include "slice.h"*/
/*^*/

#define LOOP(a) for(iy=0;iy<ayy.n;iy++){ for(ix=0;ix<axx.n;ix++){ {a} }}
#define KOOP(a) for(iy=0;iy<byy.n;iy++){ for(ix=0;ix<bxx.n;ix++){ {a} }}
#define SOOP(a) for(iy=0;iy<aly.n;iy++){ for(ix=0;ix<alx.n;ix++){ {a} }}

#define INDEX(x,a) 0.5+(x-a.o)/a.d;
#define BOUND(i,n) (i<0) ? 0 : ( (i>n-1) ? n-1 : i );
#define  KMAP(i,n) (i<n/2.) ? (i+n/2.) : (i-n/2.);
#define X2K(a,b,p) b.n=a.n+p; b.d=2.0*SF_PI/(b.n*a.d); b.o=(1==b.n)?0:-SF_PI/a.d;

static axa az,axx,ayy;
static axa    bxx,byy;
static axa    alx,aly;
static float dsmax2;

static float         **kk; /* wavenumber  */
static int            *lx;
static int            *ly;
static float complex **pk;
static float complex **wk;
static float         **wt;

/*------------------------------------------------------------*/

void ssr_init( axa az_,
	       axa ax_,
	       axa ay_,
	       axa lx_,
	       axa ly_,
	       int px,
	       int py,
	       int tx,
	       int ty,
	       float dsmax
    )
/*< initialize >*/
{
    int   ix, iy;
    float kx, ky;
    int   jx, jy;
    float xx, yy;
    int  ilx,ily;

    az  = az_;
    axx = ax_;
    ayy = ay_;
    alx = lx_;
    aly = ly_;

    /* construct K-domain axes */
    X2K(axx,bxx,px);
    X2K(ayy,byy,py);
    fft2_init(bxx.n,byy.n);

    /* precompute wavenumbers */
    kk = sf_floatalloc2(bxx.n,byy.n);  /* wavenumber */
    for (iy=0; iy<byy.n; iy++) {
	jy = KMAP(iy,byy.n);
	ky = byy.o + jy*byy.d;
	for (ix=0; ix<bxx.n; ix++) {
	    jx = KMAP(ix,bxx.n);
	    kx = bxx.o + jx*bxx.d;  
	    kk[iy][ix] = kx*kx+ky*ky;
	}
    }    

    /* precompute indices */
    lx = sf_intalloc(axx.n);
    ly = sf_intalloc(ayy.n);
    for (iy=0; iy<ayy.n; iy++) {
	yy = ayy.o + iy*ayy.d;
	ily    = INDEX( yy,aly);
	ly[iy] = BOUND(ily,aly.n); /* x-line index */
    }
    for (ix=0; ix<axx.n; ix++) {
	xx = axx.o + ix*axx.d;
	ilx    = INDEX( xx,alx);
	lx[ix] = BOUND(ilx,alx.n); /* i-line index */
    }

    /* precompute taper */
    taper2_init(ayy.n,axx.n,
		SF_MIN(ty,ayy.n-1),
		SF_MIN(tx,axx.n-1),
		true,
		true);

    /* allocate K-domain storage */
    wk = sf_complexalloc2 (bxx.n,byy.n);
    pk = sf_complexalloc2 (bxx.n,byy.n);

    /* allocate X-domain storage */
    wt = sf_floatalloc2   (axx.n,ayy.n);

    dsmax2 = dsmax*dsmax;
    dsmax2*= dsmax2;
}

/*------------------------------------------------------------*/

void ssr_close(void)
/*< free allocated storage >*/
{
    fft2_close();
    taper2_close();

    free( *kk); free( kk);
    ;           free( lx);
    ;           free( ly);

    free( *pk); free( pk);
    free( *wk); free( wk);
    free( *wt); free( wt);
}

/*------------------------------------------------------------*/

void ssr_ssf(
    float complex    w /* frequency */,
    complex float **wx /* wavefield */,
    float         **so /* slowness  */, 
    float         **ss /* slowness  */,
    int             nr /* nr. of ref slo */,
    float          *sm /* ref slo squared */
    )
/*< Wavefield extrapolation by SSF >*/
{
    float complex w2,co,cc;
    int ix,iy,jr;
    float s,d;

    w2 = w*w;

    /* w-x part 1 */
    LOOP( s = 0.5 * so[ ly[iy] ][ lx[ix] ];
	  wx[iy][ix] *= cexpf(-w*s*az.d); );

    /* FFT */
    KOOP( pk[iy][ix] = 0.; );
    LOOP( pk[iy][ix] = wx[iy][ix]; );
    fft2(false,pk);

    LOOP( wx[iy][ix] = 0;
	  wt[iy][ix] = 0; );
    for (jr=0; jr<nr; jr++) {
	/* w-k */
	co =       csqrtf(w2 * sm[jr]);
	KOOP( cc = csqrtf(w2 * sm[jr] + kk[iy][ix]);
	      wk[iy][ix] = 
	      pk[iy][ix] * cexpf((co-cc)*az.d); );

	/* IFT */
	fft2(true,wk);

	/* accumulate wavefield */
	LOOP( d = fabsf(so[ ly[iy] ][ lx[ix] ] * 
			so[ ly[iy] ][ lx[ix] ] - sm[jr]);
	      d = dsmax2/(d*d+dsmax2);
	      wx[iy][ix] += wk[iy][ix]*d;
	      wt[iy][ix] += d; );
    }
    LOOP( wx[iy][ix] /= wt[iy][ix]; );

    /* w-x part 2 */
    LOOP( s = 0.5 * ss[ ly[iy] ][ lx[ix] ];
	  wx[iy][ix] *= cexpf(-w*s*az.d); );

    taper2(wx);
}

/*------------------------------------------------------------*/

void ssr_sso(
    float complex    w /* frequency */,
    complex float **wx /* wavefield */,
    float         **so /* slowness  */, 
    float         **ss /* slowness  */,
    int             nr /* nr. of ref slo */,
    float          *sm /* ref slo squared */
    )
/*< Wavefield extrapolation by SSF w/ one reference slowness >*/
{
    float complex w2,co,cc;
    int ix,iy,jr;

    w2 = w*w;

    /* w-x part 1 */
/*    LOOP( s = 0.5 * so[ ly[iy] ][ lx[ix] ];*/
/*	  wx[iy][ix] *= cexpf(-w*s*az.d); );*/

    /* FFT */
    KOOP( pk[iy][ix] = 0.; );
    LOOP( pk[iy][ix] = wx[iy][ix]; );
    fft2(false,pk);

    jr=0;
/*    co =       csqrtf(w2 * sm[jr]);*/
    co = 0;
    KOOP( cc = csqrtf(w2 * sm[jr] + kk[iy][ix]);
	  wk[iy][ix] = 
	  pk[iy][ix] * cexpf((co-cc)*az.d); );

/*    KOOP( wk[iy][ix] = pk[iy][ix] * cexpf((co)*az.d); );*/

    /* IFT */
    fft2(true,wk);
    LOOP( wx[iy][ix] = wk[iy][ix]; );
    
    /* w-x part 2 */
/*    LOOP( s = 0.5 * ss[ ly[iy] ][ lx[ix] ];*/
/*	  wx[iy][ix] *= cexpf(-w*s*az.d); );*/

/*    taper2(wx);*/
}

/*------------------------------------------------------------*/

void ssr_phs(
    float complex    w /* frequency */,
    complex float **wx /* wavefield */,
    float         **so /* slowness  */, 
    float         **ss /* slowness  */,
    int             nr /* nr. of ref slo */,
    float          *sm /* ref slo squared */
    )
/*< Wavefield extrapolation by phase-shift >*/
{
    float complex w2,cc;
    int ix,iy,jr;
    float d;

    w2 = w*w;

    /* FFT */
    KOOP( pk[iy][ix] = 0.; );
    LOOP( pk[iy][ix] = wx[iy][ix]; );
    fft2(false,pk);

    LOOP( wx[iy][ix] = 0;
	  wt[iy][ix] = 0; );
    for (jr=0; jr<nr; jr++) {
	/* w-k */
	KOOP( cc = csqrtf(w2 * sm[jr] + kk[iy][ix]);
	      wk[iy][ix] = 
	      pk[iy][ix] * cexpf((-cc)*az.d); );

	/* IFT */
	fft2(true,wk);

	/* accumulate wavefield */
	LOOP( d = fabsf(so[ ly[iy] ][ lx[ix] ] * 
			so[ ly[iy] ][ lx[ix] ] - sm[jr]);
	      d = dsmax2/(d*d+dsmax2);
	      wx[iy][ix] += wk[iy][ix]*d;
	      wt[iy][ix] += d; );
    }
    LOOP( wx[iy][ix] /= wt[iy][ix]; );

    taper2(wx);
}
