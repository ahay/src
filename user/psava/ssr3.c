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
#include "taper3.h"

#include "weutil.h"

/*------------------------------------------------------------*/

#define LOOP(a) for(iy=0;iy<ssr->ayy.n;iy++){ for(ix=0;ix<ssr->axx.n;ix++){ {a} }}
#define KOOP(a) for(iy=0;iy<ssr->byy.n;iy++){ for(ix=0;ix<ssr->bxx.n;ix++){ {a} }}
#define SOOP(a) for(iy=0;iy<ssr->aly.n;iy++){ for(ix=0;ix<ssr->alx.n;ix++){ {a} }}

#define INDEX(x,a) 0.5+(x-a.o)/a.d;
#define BOUND(i,n) (i<0) ? 0 : ( (i>n-1) ? n-1 : i );
#define  KMAP(i,n) (i<n/2.) ? (i+n/2.) : (i-n/2.);
#define X2K(a,b,p) b.n=a.n+p; b.d=2.0*SF_PI/(b.n*a.d); b.o=(1==b.n)?0:-SF_PI/a.d;

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

    for (iy=0; iy<ssr->byy.n; iy++) {
	jy = KMAP(iy,ssr->byy.n);
	ky = ssr->byy.o + jy*ssr->byy.d;
	for (ix=0; ix<ssr->bxx.n; ix++) {
	    jx = KMAP(ix,ssr->bxx.n);
	    kx = ssr->bxx.o + jx*ssr->bxx.d;  
	    ssr->kk[iy][ix] = kx*kx+ky*ky;
	}
    }
    
    /* precompute indices */
    for (iy=0; iy<ssr->ayy.n; iy++) {
	yy = ssr->ayy.o + iy*ssr->ayy.d;
	ily    = INDEX( yy,ssr->aly);
	ssr->ly[iy] = BOUND(ily,ssr->aly.n); /* x-line index */
    }
    for (ix=0; ix<ssr->axx.n; ix++) {
	xx = ssr->axx.o + ix*ssr->axx.d;
	ilx    = INDEX( xx,ssr->alx);
	ssr->lx[ix] = BOUND(ilx,ssr->alx.n); /* i-line index */
    }

    /* initialize FFT */
    ssr->fft = ompfft2_init(ssr->bxx.n,ssr->byy.n,ssr->ompnth);

    ssr->dsmax2 = dsmax*dsmax;
    ssr->dsmax2*= ssr->dsmax2;

    return ssr;
}

/*------------------------------------------------------------*/
void ssr3_close(ssr3d ssr,
		tap3d tap)
/*< free allocated storage >*/
{
    ompfft2_close(ssr->fft);
    taper2d_close(tap);

    free( *ssr->kk); free( ssr->kk);
    ;           free( ssr->lx);
    ;           free( ssr->ly);

    free(**ssr->pk); free( *ssr->pk); free( ssr->pk);
    free(**ssr->wk); free( *ssr->wk); free( ssr->wk);
    free(**ssr->wt); free( *ssr->wt); free( ssr->wt);
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
    ssr3d          ssr,
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
    LOOP( s = 0.5 * so[ ssr->ly[iy] ][ ssr->lx[ix] ];
	  wx[iy][ix] *= cexpf(-w*s*ssr->az.d); );
#else
    w2 = sf_cmul(w,w);
    /* w-x part 1 */
    LOOP( s = 0.5 * so[ ssr->ly[iy] ][ ssr->lx[ix] ];
	  wx[iy][ix] = sf_cmul( wx[iy][ix],cexpf(sf_crmul(w,-s*ssr->az.d))); );
#endif

    /* FFT */
    KOOP( ssr->pk[ompith][iy][ix] = sf_cmplx(0.,0.); );
    LOOP( ssr->pk[ompith][iy][ix] = wx[iy][ix]; );
    ompfft2(false,(kiss_fft_cpx**) ssr->pk[ompith],ompith,ssr->fft);

    LOOP( wx[iy][ix] = sf_cmplx(0.,0.);
	  ssr->wt[ompith][iy][ix] = 0; );
    for (jr=0; jr<nr; jr++) {
	/* w-k */
#ifdef SF_HAS_COMPLEX_H
	co =       csqrtf(w2 * sm[jr]);
	KOOP( cc = csqrtf(w2 * sm[jr] + ssr->kk[iy][ix]);
	      ssr->wk[ompith][iy][ix] = 
	      ssr->pk[ompith][iy][ix] * cexpf((co-cc)*ssr->az.d); 
	    );
#else
	co =       csqrtf(sf_crmul(w2,sm[jr]));
	KOOP( cc = csqrtf(sf_cadd(sf_crmul(w2,sm[jr]),
				  sf_cmplx(ssr->kk[iy][ix],0.)));
	      ssr->wk[ompith][iy][ix] = 
	      sf_cmul(pk[ompith][iy][ix],cexpf(sf_crmul(sf_csub(co,cc),ssr->az.d))); 
	    );
#endif
	
	/* IFT */
	ompfft2(true,(kiss_fft_cpx**) ssr->wk[ompith],ompith,ssr->fft);
	
	/* accumulate wavefield */
#ifdef SF_HAS_COMPLEX_H
	LOOP( d = fabsf(so[ ssr->ly[iy] ][ ssr->lx[ix] ] * 
			so[ ssr->ly[iy] ][ ssr->lx[ix] ] - sm[jr]);
	      d = ssr->dsmax2/(d*d+ssr->dsmax2);
	      wx[iy][ix] += ssr->wk[ompith][iy][ix]*d;
	      ssr->wt[ompith][iy][ix] += d; 
	    );
#else
	LOOP( d = fabsf(so[ ssr->ly[iy] ][ ssr->lx[ix] ] * 
			so[ ssr->ly[iy] ][ ssr->lx[ix] ] - sm[jr]);
	      d = ssr->dsmax2/(d*d+ssr->dsmax2);
	      wx[iy][ix] = sf_cadd(wx[iy][ix],sf_crmul(ssr->wk[ompith][iy][ix],d));
	      ssr->wt[ompith][iy][ix] += d; 
	    );
#endif
    }
#ifdef SF_HAS_COMPLEX_H
    LOOP( wx[iy][ix] /= ssr->wt[ompith][iy][ix]; );
    /* w-x part 2 */
    LOOP( s = 0.5 * ss[ ssr->ly[iy] ][ ssr->lx[ix] ];
	  wx[iy][ix] *= cexpf(-w*s*ssr->az.d); );
#else
    LOOP( wx[iy][ix] = sf_crmul(wx[iy][ix],1./ssr->wt[ompith][iy][ix]); );
    /* w-x part 2 */
    LOOP( s = 0.5 * ss[ ssr->ly[iy] ][ ssr->lx[ix] ];
	  wx[iy][ix] = sf_cmul(wx[iy][ix],cexpf(sf_crmul(w,-s*ssr->az.d))); );
#endif

    taper2d(wx,tap);
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
    ssr3d          ssr,
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
    LOOP( s = 0.5 * so[ ssr->ly[iy] ][ ssr->lx[ix] ];
	  wx[iy][ix] *= cexpf(-w*s*ssr->az.d); );
#else
    w2 = sf_cmul(w,w);
    /* w-x part 1 */
    LOOP( s = 0.5 * so[ ssr->ly[iy] ][ ssr->lx[ix] ];
	  wx[iy][ix] = sf_cmul(wx[iy][ix],cexpf(sf_crmul(w,-s*ssr->az.d))); );
#endif

    /* FFT */
    KOOP( ssr->pk[ompith][iy][ix] = sf_cmplx(0.,0.); );
    LOOP( ssr->pk[ompith][iy][ix] = wx[iy][ix]; );
    ompfft2(false,(kiss_fft_cpx**) ssr->pk[ompith],ompith,ssr->fft);

    jr=0;
#ifdef SF_HAS_COMPLEX_H
    co =       csqrtf(w2 * sm[jr]);
    KOOP( cc = csqrtf(w2 * sm[jr] + ssr->kk[iy][ix]);
	  ssr->wk[ompith][iy][ix] = 
	  ssr->pk[ompith][iy][ix] * cexpf((co-cc)*ssr->az.d); 
	);

    KOOP( ssr->wk[ompith][iy][ix] = 
	  ssr->pk[ompith][iy][ix] * cexpf((co)*ssr->az.d); );
#else
    co =       csqrtf(sf_crmul(w2,sm[jr]));
    KOOP( cc = csqrtf(sf_cadd(sf_crmul(w2,sm[jr]),
			      sf_cmplx(ssr->kk[iy][ix],0.)));
	  ssr->wk[ompith][iy][ix] = sf_cmul(ssr->pk[ompith][iy][ix],
			       cexpf(sf_crmul(sf_csub(co,cc),ssr->az.d))); 
	);

    KOOP( ssr->wk[ompith][iy][ix] = sf_cmul(ssr->pk[ompith][iy][ix],cexpf(sf_crmul(co,ssr->az.d))); );
#endif

    /* IFT */
    ompfft2(true,(kiss_fft_cpx**) ssr->wk[ompith],ompith,ssr->fft);
    LOOP( wx[iy][ix] = ssr->wk[ompith][iy][ix]; );
    
    /* w-x part 2 */
#ifdef SF_HAS_COMPLEX_H
    LOOP( s = 0.5 * ss[ ssr->ly[iy] ][ ssr->lx[ix] ];
	  wx[iy][ix] *= cexpf(-w*s*ssr->az.d); );
#else
    LOOP( s = 0.5 * ss[ ssr->ly[iy] ][ ssr->lx[ix] ];
	  wx[iy][ix] = sf_cmul(wx[iy][ix],cexpf(sf_crmul(w,-s*ssr->az.d))); );
#endif

    taper2d(wx,tap);
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
    ssr3d          ssr,
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
    KOOP( ssr->pk[ompith][iy][ix] = sf_cmplx(0.,0.); );
    LOOP( ssr->pk[ompith][iy][ix] = wx[iy][ix]; );
    ompfft2(false,(kiss_fft_cpx**) ssr->pk[ompith],ompith,ssr->fft);

    LOOP( wx[iy][ix] = sf_cmplx(0.,0.);
	  ssr->wt[ompith][iy][ix] = 0; );
    for (jr=0; jr<nr; jr++) {
	/* w-k */
#ifdef SF_HAS_COMPLEX_H
	KOOP( cc = csqrtf(w2 * sm[jr] + ssr->kk[iy][ix]);
	      ssr->wk[ompith][iy][ix] = 
	      ssr->pk[ompith][iy][ix] * cexpf((-cc)*ssr->az.d); 
	    );
#else
	KOOP( cc = csqrtf(sf_cadd(sf_crmul(w2,sm[jr]),
				  sf_cmplx(ssr->kk[iy][ix],0.)));
	      ssr->wk[ompith][iy][ix] = 
	      sf_cmul(ssr->pk[ompith][iy][ix],cexpf(sf_crmul(cc,-ssr->az.d))); 
	    );
#endif

	/* IFT */
	ompfft2(true,(kiss_fft_cpx**) ssr->wk[ompith],ompith,ssr->fft);

	/* accumulate wavefield */
#ifdef SF_HAS_COMPLEX_H
	LOOP( d = fabsf(so[ ssr->ly[iy] ][ ssr->lx[ix] ] * 
			so[ ssr->ly[iy] ][ ssr->lx[ix] ] - sm[jr]);
	      d = ssr->dsmax2/(d*d+ssr->dsmax2);
	      wx[iy][ix] += ssr->wk[ompith][iy][ix]*d;
	      ssr->wt[ompith][iy][ix] += d; );
#else
	LOOP( d = fabsf(so[ ssr->ly[iy] ][ ssr->lx[ix] ] * 
			so[ ssr->ly[iy] ][ ssr->lx[ix] ] - sm[jr]);
	      d = ssr->dsmax2/(d*d+ssr->dsmax2);
	      wx[iy][ix] = sf_cadd(wx[iy][ix],sf_crmul(ssr->wk[ompith][iy][ix],d));
	      ssr->wt[ompith][iy][ix] += d; );
#endif
    }
#ifdef SF_HAS_COMPLEX_H
    LOOP( wx[iy][ix] /= ssr->wt[ompith][iy][ix]; );
#else
    LOOP( wx[iy][ix] = sf_crmul(wx[iy][ix],1./ssr->wt[ompith][iy][ix]); );
#endif

    taper2d(wx,tap);
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
    ssr3d          ssr,
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
    KOOP( ssr->pk[ompith][iy][ix] = sf_cmplx(0.,0.); );
    LOOP( ssr->pk[ompith][iy][ix] = wx[iy][ix]; );
    ompfft2(false,(kiss_fft_cpx**) ssr->pk[ompith],ompith,ssr->fft);

    jr=0;
#ifdef SF_HAS_COMPLEX_H
    KOOP(
	cc = csqrtf(w2 * sm[jr] + ssr->kk[iy][ix]);
	
	ssr->wk[ompith][iy][ix] = 
	ssr->pk[ompith][iy][ix] * cexpf(-cc*ssr->az.d); 
	);
#else
    KOOP(
	cc = csqrtf(sf_cadd(sf_crmul(w2,sm[jr]),
			    sf_cmplx(ssr->kk[iy][ix],0.)));
	
	ssr->wk[ompith][iy][ix] = 
	sf_cmul(ssr->pk[ompith][iy][ix],cexpf(sf_crmul(cc,-ssr->az.d))); 
	);
#endif
    
    /* IFT */
    ompfft2( true,(kiss_fft_cpx**) ssr->wk[ompith],ompith,ssr->fft);    
    LOOP( wx[iy][ix] = ssr->wk[ompith][iy][ix]; );
    
    taper2d(wx,tap);
    
}

