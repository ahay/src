/* 3-D SSR */

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
/*^*/

#include "fft2.h"
#include "taper.h"

#define LOOP(a) for(iy=0;iy<ayy.n;iy++){ for(ix=0;ix<axx.n;ix++){ {a} }}
#define KOOP(a) for(iy=0;iy<byy.n;iy++){ for(ix=0;ix<bxx.n;ix++){ {a} }}
#define SOOP(a) for(iy=0;iy<aly.n;iy++){ for(ix=0;ix<alx.n;ix++){ {a} }}

#define INDEX(x,a) 0.5+(x-a.o)/a.d;
#define BOUND(i,n) (i<0) ? 0 : ( (i>n-1) ? n-1 : i );
#define  KMAP(i,n) (i<n/2.) ? (i+n/2.) : (i-n/2.);
#define X2K(a,b,p) b.n=a.n+p; b.d=2.0*SF_PI/(b.n*a.d); b.o=(1==b.n)?0:-SF_PI/a.d;

static sf_axa az,axx,ayy;
static sf_axa    bxx,byy;
static sf_axa    alx,aly;
static float dsmax2;

static float         **kk; /* wavenumber  */
static int            *lx;
static int            *ly;
static sf_complex **pk;
static sf_complex **wk;
static float         **wt;

/*------------------------------------------------------------*/
void ssr_init( sf_axis az_,
	       sf_axis ax_,
	       sf_axis ay_,
	       sf_axis lx_,
	       sf_axis ly_,
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

    az  = sf_nod(az_);
    axx = sf_nod(ax_);
    ayy = sf_nod(ay_);
    alx = sf_nod(lx_);
    aly = sf_nod(ly_);

    /* construct K-domain axes */
    X2K(axx,bxx,px); sf_raxa(sf_maxa(bxx.n,bxx.o,bxx.d));
    X2K(ayy,byy,py); sf_raxa(sf_maxa(byy.n,byy.o,byy.d));
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

/*	    sf_warning("kk[%d][%d]=%g",iy,ix,kk[iy][ix]);*/
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
    sf_complex    w /* frequency */,
    sf_complex **wx /* wavefield */,
    float         **so /* slowness  */, 
    float         **ss /* slowness  */,
    int             nr /* nr. of ref slo */,
    float          *sm /* ref slo squared */
    )
/*< Wavefield extrapolation by SSF >*/
{
    sf_complex w2,co,cc;
    int ix,iy,jr;
    float s,d;

#ifdef SF_HAS_COMPLEX_H
    w2 = w*w;
    /* w-x part 1 */
    LOOP( s = 0.5 * so[ ly[iy] ][ lx[ix] ];
	  wx[iy][ix] *= cexpf(-w*s*az.d); );
#else
    w2 = sf_cmul(w,w);
    /* w-x part 1 */
    LOOP( s = 0.5 * so[ ly[iy] ][ lx[ix] ];
	  wx[iy][ix] = sf_cmul( wx[iy][ix],cexpf(sf_crmul(w,-s*az.d))); );
#endif

    /* FFT */
    KOOP( pk[iy][ix] = sf_cmplx(0.,0.); );
    LOOP( pk[iy][ix] = wx[iy][ix]; );
    fft2(false,(kiss_fft_cpx**) pk);

    LOOP( wx[iy][ix] = sf_cmplx(0.,0.);
	  wt[iy][ix] = 0; );
    for (jr=0; jr<nr; jr++) {
	/* w-k */
#ifdef SF_HAS_COMPLEX_H
	co =       csqrtf(w2 * sm[jr]);
	KOOP( cc = csqrtf(w2 * sm[jr] + kk[iy][ix]);
	      wk[iy][ix] = 
	      pk[iy][ix] * cexpf((co-cc)*az.d); 
	    );
#else
	co =       csqrtf(sf_crmul(w2,sm[jr]));
	KOOP( cc = csqrtf(sf_cadd(sf_crmul(w2,sm[jr]),
				  sf_cmplx(kk[iy][ix],0.)));
	      wk[iy][ix] = 
	      sf_cmul(pk[iy][ix],cexpf(sf_crmul(sf_csub(co,cc),az.d))); 
	    );
#endif
	
	/* IFT */
	fft2(true,(kiss_fft_cpx**) wk);
	
	/* accumulate wavefield */
#ifdef SF_HAS_COMPLEX_H
	LOOP( d = fabsf(so[ ly[iy] ][ lx[ix] ] * 
			so[ ly[iy] ][ lx[ix] ] - sm[jr]);
	      d = dsmax2/(d*d+dsmax2);
	      wx[iy][ix] += wk[iy][ix]*d;
	      wt[iy][ix] += d; 
	    );
#else
	LOOP( d = fabsf(so[ ly[iy] ][ lx[ix] ] * 
			so[ ly[iy] ][ lx[ix] ] - sm[jr]);
	      d = dsmax2/(d*d+dsmax2);
	      wx[iy][ix] = sf_cadd(wx[iy][ix],sf_crmul(wk[iy][ix],d));
	      wt[iy][ix] += d; 
	    );
#endif
    }
#ifdef SF_HAS_COMPLEX_H
    LOOP( wx[iy][ix] /= wt[iy][ix]; );
    /* w-x part 2 */
    LOOP( s = 0.5 * ss[ ly[iy] ][ lx[ix] ];
	  wx[iy][ix] *= cexpf(-w*s*az.d); );
#else
    LOOP( wx[iy][ix] = sf_crmul(wx[iy][ix],1./wt[iy][ix]); );
    /* w-x part 2 */
    LOOP( s = 0.5 * ss[ ly[iy] ][ lx[ix] ];
	  wx[iy][ix] = sf_cmul(wx[iy][ix],cexpf(sf_crmul(w,-s*az.d))); );
#endif

    taper2(wx);
}

/*------------------------------------------------------------*/
void ssr_sso(
    sf_complex    w /* frequency */,
    sf_complex **wx /* wavefield */,
    float      **so /* slowness  */, 
    float      **ss /* slowness  */,
    int          nr /* nr. of ref slo */,
    float       *sm /* ref slo squared */
    )
/*< Wavefield extrapolation by SSF w/ one reference slowness >*/
{
    sf_complex w2,co,cc;
    int ix,iy,jr;
    float s;

#ifdef SF_HAS_COMPLEX_H
    w2 = w*w;
    /* w-x part 1 */
    LOOP( s = 0.5 * so[ ly[iy] ][ lx[ix] ];
	  wx[iy][ix] *= cexpf(-w*s*az.d); );
#else
    w2 = sf_cmul(w,w);
    /* w-x part 1 */
    LOOP( s = 0.5 * so[ ly[iy] ][ lx[ix] ];
	  wx[iy][ix] = sf_cmul(wx[iy][ix],cexpf(sf_crmul(w,-s*az.d))); );
#endif

    /* FFT */
    KOOP( pk[iy][ix] = sf_cmplx(0.,0.); );
    LOOP( pk[iy][ix] = wx[iy][ix]; );
    fft2(false,(kiss_fft_cpx**) pk);

    jr=0;
#ifdef SF_HAS_COMPLEX_H
    co =       csqrtf(w2 * sm[jr]);
    KOOP( cc = csqrtf(w2 * sm[jr] + kk[iy][ix]);
	  wk[iy][ix] = 
	  pk[iy][ix] * cexpf((co-cc)*az.d); 
	);

    KOOP( wk[iy][ix] = pk[iy][ix] * cexpf((co)*az.d); );
#else
    co =       csqrtf(sf_crmul(w2,sm[jr]));
    KOOP( cc = csqrtf(sf_cadd(sf_crmul(w2,sm[jr]),
			      sf_cmplx(kk[iy][ix],0.)));
	  wk[iy][ix] = sf_cmul(pk[iy][ix],
			       cexpf(sf_crmul(sf_csub(co,cc),az.d))); 
	);

    KOOP( wk[iy][ix] = sf_cmul(pk[iy][ix],cexpf(sf_crmul(co,az.d))); );
#endif

    /* IFT */
    fft2(true,(kiss_fft_cpx**) wk);
    LOOP( wx[iy][ix] = wk[iy][ix]; );
    
    /* w-x part 2 */
#ifdef SF_HAS_COMPLEX_H
    LOOP( s = 0.5 * ss[ ly[iy] ][ lx[ix] ];
	  wx[iy][ix] *= cexpf(-w*s*az.d); );
#else
    LOOP( s = 0.5 * ss[ ly[iy] ][ lx[ix] ];
	  wx[iy][ix] = sf_cmul(wx[iy][ix],cexpf(sf_crmul(w,-s*az.d))); );
#endif

    taper2(wx);
}

/*------------------------------------------------------------*/
void ssr_phs(
    sf_complex    w /* frequency */,
    sf_complex **wx /* wavefield */,
    float      **so /* slowness  */, 
    float      **ss /* slowness  */,
    int          nr /* nr. of ref slo */,
    float       *sm /* ref slo squared */
    )
/*< Wavefield extrapolation by phase-shift >*/
{
    sf_complex w2,cc;
    int ix,iy,jr;
    float d;

#ifdef SF_HAS_COMPLEX_H
    w2 = w*w;
#else
    w2 = sf_cmul(w,w);
#endif

    /* FFT */
    KOOP( pk[iy][ix] = sf_cmplx(0.,0.); );
    LOOP( pk[iy][ix] = wx[iy][ix]; );
    fft2(false,(kiss_fft_cpx**) pk);

    LOOP( wx[iy][ix] = sf_cmplx(0.,0.);
	  wt[iy][ix] = 0; );
    for (jr=0; jr<nr; jr++) {
	/* w-k */
#ifdef SF_HAS_COMPLEX_H
	KOOP( cc = csqrtf(w2 * sm[jr] + kk[iy][ix]);
	      wk[iy][ix] = 
	      pk[iy][ix] * cexpf((-cc)*az.d); 
	    );
#else
	KOOP( cc = csqrtf(sf_cadd(sf_crmul(w2,sm[jr]),
				  sf_cmplx(kk[iy][ix],0.)));
	      wk[iy][ix] = 
	      sf_cmul(pk[iy][ix],cexpf(sf_crmul(cc,-az.d))); 
	    );
#endif

	/* IFT */
	fft2(true,(kiss_fft_cpx**) wk);

	/* accumulate wavefield */
#ifdef SF_HAS_COMPLEX_H
	LOOP( d = fabsf(so[ ly[iy] ][ lx[ix] ] * 
			so[ ly[iy] ][ lx[ix] ] - sm[jr]);
	      d = dsmax2/(d*d+dsmax2);
	      wx[iy][ix] += wk[iy][ix]*d;
	      wt[iy][ix] += d; );
#else
	LOOP( d = fabsf(so[ ly[iy] ][ lx[ix] ] * 
			so[ ly[iy] ][ lx[ix] ] - sm[jr]);
	      d = dsmax2/(d*d+dsmax2);
	      wx[iy][ix] = sf_cadd(wx[iy][ix],sf_crmul(wk[iy][ix],d));
	      wt[iy][ix] += d; );
#endif
    }
#ifdef SF_HAS_COMPLEX_H
    LOOP( wx[iy][ix] /= wt[iy][ix]; );
#else
    LOOP( wx[iy][ix] = sf_crmul(wx[iy][ix],1./wt[iy][ix]); );
#endif

    taper2(wx);
}

/*------------------------------------------------------------*/
void ssr_pho(
    sf_complex    w /* frequency */,
    sf_complex **wx /* wavefield */,
    float      **so /* slowness  */, 
    float      **ss /* slowness  */,
    int          nr /* nr. of ref slo */,
    float       *sm /* ref slo squared */
    )
/*< Wavefield extrapolation by phase-shift w/ one reference >*/
{
    sf_complex  w2,cc;
    int ix,iy,jr;

    w = sf_cmplx(1e-3,cimagf(w));

#ifdef SF_HAS_COMPLEX_H
    w2 = w*w;
#else
    w2 = sf_cmul(w,w);
#endif

    /* FFT */
    KOOP( pk[iy][ix] = sf_cmplx(0.,0.); );
    LOOP( pk[iy][ix] = wx[iy][ix]; );
    fft2(false,(kiss_fft_cpx**) pk);

    jr=0;
#ifdef SF_HAS_COMPLEX_H
    KOOP(
	cc = csqrtf(w2 * sm[jr] + kk[iy][ix]);

	wk[iy][ix] = 
	pk[iy][ix] * cexpf(-cc*az.d); 
	);
#else
     KOOP(
	cc = csqrtf(sf_cadd(sf_crmul(w2,sm[jr]),
			    sf_cmplx(kk[iy][ix],0.)));

	wk[iy][ix] = 
	sf_cmul(pk[iy][ix],cexpf(sf_crmul(cc,-az.d))); 
	);
#endif

    /* IFT */
    fft2( true,(kiss_fft_cpx**) wk);
    LOOP( wx[iy][ix] = wk[iy][ix]; );
    
    taper2(wx);
}

