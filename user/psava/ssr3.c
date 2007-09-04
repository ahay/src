/* 3-D SSR */

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

#include "ompfft2.h"
#include "taper.h"

#include "weutil.h"

/*------------------------------------------------------------*/

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
static sf_complex   ***pk;
static sf_complex   ***wk;
static float        ***wt;

/*------------------------------------------------------------*/
ssr3d ssr3_init( sf_axis az_,
		 sf_axis ax_,
		 sf_axis ay_,
		 sf_axis lx_,
		 sf_axis ly_,
		 int px,
		 int py,
		 int tx,
		 int ty,
		 float dsmax,
		 int ompnth
    )
/*< initialize >*/
{
    int   ix, iy;
    float kx, ky;
    int   jx, jy;
    float xx, yy;
    int  ilx,ily;


    /*------------------------------------------------------------*/

    ssr3d ssr;
    ssr = (ssr3d) sf_alloc(1,sizeof(*ssr));

    ssr->ompnth = ompnth;

    ssr->az  = sf_nod(az_);
    ssr->axx = sf_nod(ax_);
    ssr->ayy = sf_nod(ay_);
    ssr->alx = sf_nod(lx_);
    ssr->aly = sf_nod(ly_);

    X2K(ssr->axx,ssr->bxx,px);
    X2K(ssr->ayy,ssr->byy,py);

    ssr->kk = sf_floatalloc2(ssr->bxx.n,ssr->byy.n);

    ssr->lx = sf_intalloc(ssr->axx.n);
    ssr->ly = sf_intalloc(ssr->ayy.n);

    /* allocate K-domain storage */
    ssr->wk = sf_complexalloc3 (ssr->bxx.n,ssr->byy.n,ssr->ompnth);
    ssr->pk = sf_complexalloc3 (ssr->bxx.n,ssr->byy.n,ssr->ompnth);

    /* allocate X-domain storage */
    ssr->wt = sf_floatalloc3   (ssr->axx.n,ssr->ayy.n,ssr->ompnth);

    /*------------------------------------------------------------*/


    az  = sf_nod(az_);
    axx = sf_nod(ax_);
    ayy = sf_nod(ay_);
    alx = sf_nod(lx_);
    aly = sf_nod(ly_);

    /* construct K-domain axes */
    X2K(axx,bxx,px); sf_raxa(sf_maxa(bxx.n,bxx.o,bxx.d));
    X2K(ayy,byy,py); sf_raxa(sf_maxa(byy.n,byy.o,byy.d));

    /* allocate K-domain storage */
    wk = sf_complexalloc3 (bxx.n,byy.n,ompnth);
    pk = sf_complexalloc3 (bxx.n,byy.n,ompnth);

    /* allocate X-domain storage */
    wt = sf_floatalloc3   (axx.n,ayy.n,ompnth);


    ompfft2_init(bxx.n,byy.n,ompnth);

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

    for (iy=0; iy<byy.n; iy++) {
	jy = KMAP(iy,ssr->byy.n);
	ky = ssr->byy.o + jy*ssr->byy.d;
	for (ix=0; ix<ssr->bxx.n; ix++) {
	    jx = KMAP(ix,ssr->bxx.n);
	    kx = ssr->bxx.o + jx*ssr->bxx.d;  
	    ssr->kk[iy][ix] = kx*kx+ky*ky;
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

    for (iy=0; iy<ssr->ayy.n; iy++) {
	yy = ssr->ayy.o + iy*ayy.d;
	ily    = INDEX( yy,ssr->aly);
	ssr->ly[iy] = BOUND(ily,ssr->aly.n); /* x-line index */
    }
    for (ix=0; ix<ssr->axx.n; ix++) {
	xx = ssr->axx.o + ix*ssr->axx.d;
	ilx    = INDEX( xx,ssr->alx);
	ssr->lx[ix] = BOUND(ilx,ssr->alx.n); /* i-line index */
    }

    /* initialize FFT */
    ompfft2_init(ssr->bxx.n,ssr->byy.n,ssr->ompnth);

    /* precompute taper */
/*    taper2_init(ssr->ayy.n,ssr->axx.n,*/
/*		SF_MIN(ty,ssr->ayy.n-1),*/
/*		SF_MIN(tx,ssr->axx.n-1),*/
/*		true,*/
/*		true);*/

    dsmax2 = dsmax*dsmax;
    dsmax2*= dsmax2;

    ssr->dsmax2 = dsmax*dsmax;
    ssr->dsmax2*= ssr->dsmax2;

    return ssr;
}

/*------------------------------------------------------------*/
void ssr3_close(tap3d tap)
/*< free allocated storage >*/
{
    ompfft2_close();
    taper2_close(tap);

    free( *kk); free( kk);
    ;           free( lx);
    ;           free( ly);

    free(**pk); free( *pk); free( pk);
    free(**wk); free( *wk); free( wk);
    free(**wt); free( *wt); free( wt);
}

/*------------------------------------------------------------*/
void ssr3_ssf(
    sf_complex       w /* frequency */,
    sf_complex    **wx /* wavefield */,
    float         **so /* slowness  */, 
    float         **ss /* slowness  */,
    int             nr /* nr. of ref slo */,
    float          *sm /* ref slo squared */,
    int         ompith,
    tap3d          tap
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
    KOOP( pk[ompith][iy][ix] = sf_cmplx(0.,0.); );
    LOOP( pk[ompith][iy][ix] = wx[iy][ix]; );
    ompfft2(false,(kiss_fft_cpx**) pk[ompith],ompith);

    LOOP( wx[iy][ix] = sf_cmplx(0.,0.);
	  wt[ompith][iy][ix] = 0; );
    for (jr=0; jr<nr; jr++) {
	/* w-k */
#ifdef SF_HAS_COMPLEX_H
	co =       csqrtf(w2 * sm[jr]);
	KOOP( cc = csqrtf(w2 * sm[jr] + kk[iy][ix]);
	      wk[ompith][iy][ix] = 
	      pk[ompith][iy][ix] * cexpf((co-cc)*az.d); 
	    );
#else
	co =       csqrtf(sf_crmul(w2,sm[jr]));
	KOOP( cc = csqrtf(sf_cadd(sf_crmul(w2,sm[jr]),
				  sf_cmplx(kk[iy][ix],0.)));
	      wk[ompith][iy][ix] = 
	      sf_cmul(pk[ompith][iy][ix],cexpf(sf_crmul(sf_csub(co,cc),az.d))); 
	    );
#endif
	
	/* IFT */
	ompfft2(true,(kiss_fft_cpx**) wk[ompith],ompith);
	
	/* accumulate wavefield */
#ifdef SF_HAS_COMPLEX_H
	LOOP( d = fabsf(so[ ly[iy] ][ lx[ix] ] * 
			so[ ly[iy] ][ lx[ix] ] - sm[jr]);
	      d = dsmax2/(d*d+dsmax2);
	      wx[iy][ix] += wk[ompith][iy][ix]*d;
	      wt[ompith][iy][ix] += d; 
	    );
#else
	LOOP( d = fabsf(so[ ly[iy] ][ lx[ix] ] * 
			so[ ly[iy] ][ lx[ix] ] - sm[jr]);
	      d = dsmax2/(d*d+dsmax2);
	      wx[iy][ix] = sf_cadd(wx[iy][ix],sf_crmul(wk[ompith][iy][ix],d));
	      wt[ompith][iy][ix] += d; 
	    );
#endif
    }
#ifdef SF_HAS_COMPLEX_H
    LOOP( wx[iy][ix] /= wt[ompith][iy][ix]; );
    /* w-x part 2 */
    LOOP( s = 0.5 * ss[ ly[iy] ][ lx[ix] ];
	  wx[iy][ix] *= cexpf(-w*s*az.d); );
#else
    LOOP( wx[iy][ix] = sf_crmul(wx[iy][ix],1./wt[ompith][iy][ix]); );
    /* w-x part 2 */
    LOOP( s = 0.5 * ss[ ly[iy] ][ lx[ix] ];
	  wx[iy][ix] = sf_cmul(wx[iy][ix],cexpf(sf_crmul(w,-s*az.d))); );
#endif

    taper2(wx,tap);
}

/*------------------------------------------------------------*/
void ssr3_sso(
    sf_complex    w /* frequency */,
    sf_complex **wx /* wavefield */,
    float      **so /* slowness  */, 
    float      **ss /* slowness  */,
    int          nr /* nr. of ref slo */,
    float       *sm /* ref slo squared */,
    int         ompith,
    tap3d          tap
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
    KOOP( pk[ompith][iy][ix] = sf_cmplx(0.,0.); );
    LOOP( pk[ompith][iy][ix] = wx[iy][ix]; );
    ompfft2(false,(kiss_fft_cpx**) pk[ompith],ompith);

    jr=0;
#ifdef SF_HAS_COMPLEX_H
    co =       csqrtf(w2 * sm[jr]);
    KOOP( cc = csqrtf(w2 * sm[jr] + kk[iy][ix]);
	  wk[ompith][iy][ix] = 
	  pk[ompith][iy][ix] * cexpf((co-cc)*az.d); 
	);

    KOOP( wk[ompith][iy][ix] = pk[ompith][iy][ix] * cexpf((co)*az.d); );
#else
    co =       csqrtf(sf_crmul(w2,sm[jr]));
    KOOP( cc = csqrtf(sf_cadd(sf_crmul(w2,sm[jr]),
			      sf_cmplx(kk[iy][ix],0.)));
	  wk[ompith][iy][ix] = sf_cmul(pk[ompith][iy][ix],
			       cexpf(sf_crmul(sf_csub(co,cc),az.d))); 
	);

    KOOP( wk[ompith][iy][ix] = sf_cmul(pk[ompith][iy][ix],cexpf(sf_crmul(co,az.d))); );
#endif

    /* IFT */
    ompfft2(true,(kiss_fft_cpx**) wk[ompith],ompith);
    LOOP( wx[iy][ix] = wk[ompith][iy][ix]; );
    
    /* w-x part 2 */
#ifdef SF_HAS_COMPLEX_H
    LOOP( s = 0.5 * ss[ ly[iy] ][ lx[ix] ];
	  wx[iy][ix] *= cexpf(-w*s*az.d); );
#else
    LOOP( s = 0.5 * ss[ ly[iy] ][ lx[ix] ];
	  wx[iy][ix] = sf_cmul(wx[iy][ix],cexpf(sf_crmul(w,-s*az.d))); );
#endif

    taper2(wx,tap);
}

/*------------------------------------------------------------*/
void ssr3_phs(
    sf_complex    w /* frequency */,
    sf_complex **wx /* wavefield */,
    float      **so /* slowness  */, 
    float      **ss /* slowness  */,
    int          nr /* nr. of ref slo */,
    float       *sm /* ref slo squared */,
    int         ompith,
    tap3d          tap
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
    KOOP( pk[ompith][iy][ix] = sf_cmplx(0.,0.); );
    LOOP( pk[ompith][iy][ix] = wx[iy][ix]; );
    ompfft2(false,(kiss_fft_cpx**) pk[ompith],ompith);

    LOOP( wx[iy][ix] = sf_cmplx(0.,0.);
	  wt[ompith][iy][ix] = 0; );
    for (jr=0; jr<nr; jr++) {
	/* w-k */
#ifdef SF_HAS_COMPLEX_H
	KOOP( cc = csqrtf(w2 * sm[jr] + kk[iy][ix]);
	      wk[ompith][iy][ix] = 
	      pk[ompith][iy][ix] * cexpf((-cc)*az.d); 
	    );
#else
	KOOP( cc = csqrtf(sf_cadd(sf_crmul(w2,sm[jr]),
				  sf_cmplx(kk[iy][ix],0.)));
	      wk[ompith][iy][ix] = 
	      sf_cmul(pk[ompith][iy][ix],cexpf(sf_crmul(cc,-az.d))); 
	    );
#endif

	/* IFT */
	ompfft2(true,(kiss_fft_cpx**) wk[ompith],ompith);

	/* accumulate wavefield */
#ifdef SF_HAS_COMPLEX_H
	LOOP( d = fabsf(so[ ly[iy] ][ lx[ix] ] * 
			so[ ly[iy] ][ lx[ix] ] - sm[jr]);
	      d = dsmax2/(d*d+dsmax2);
	      wx[iy][ix] += wk[ompith][iy][ix]*d;
	      wt[ompith][iy][ix] += d; );
#else
	LOOP( d = fabsf(so[ ly[iy] ][ lx[ix] ] * 
			so[ ly[iy] ][ lx[ix] ] - sm[jr]);
	      d = dsmax2/(d*d+dsmax2);
	      wx[iy][ix] = sf_cadd(wx[iy][ix],sf_crmul(wk[ompith][iy][ix],d));
	      wt[ompith][iy][ix] += d; );
#endif
    }
#ifdef SF_HAS_COMPLEX_H
    LOOP( wx[iy][ix] /= wt[ompith][iy][ix]; );
#else
    LOOP( wx[iy][ix] = sf_crmul(wx[iy][ix],1./wt[ompith][iy][ix]); );
#endif

    taper2(wx,tap);
}

/*------------------------------------------------------------*/
void ssr3_pho(
    sf_complex    w /* frequency */,
    sf_complex **wx /* wavefield */,
    float      **so /* slowness  */, 
    float      **ss /* slowness  */,
    int          nr /* nr. of ref slo */,
    float       *sm /* ref slo squared */,
    int         ompith,
    tap3d          tap
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
    KOOP( pk[ompith][iy][ix] = sf_cmplx(0.,0.); );
    LOOP( pk[ompith][iy][ix] = wx[iy][ix]; );
    ompfft2(false,(kiss_fft_cpx**) pk[ompith],ompith);

    jr=0;
#ifdef SF_HAS_COMPLEX_H
    KOOP(
	cc = csqrtf(w2 * sm[jr] + kk[iy][ix]);
	
	wk[ompith][iy][ix] = 
	pk[ompith][iy][ix] * cexpf(-cc*az.d); 
	);
#else
    KOOP(
	cc = csqrtf(sf_cadd(sf_crmul(w2,sm[jr]),
			    sf_cmplx(kk[iy][ix],0.)));
	
	wk[ompith][iy][ix] = 
	sf_cmul(pk[ompith][iy][ix],cexpf(sf_crmul(cc,-az.d))); 
	);
#endif
    
    /* IFT */
    ompfft2( true,(kiss_fft_cpx**) wk[ompith],ompith);    
    LOOP( wx[iy][ix] = wk[ompith][iy][ix]; );
    
    taper2(wx,tap);
    
}

