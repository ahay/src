/* 3-D SSR using extended SSF */

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

#define LOOP(a) for(iy=0;iy<cub->amy.n;iy++){ \
	        for(ix=0;ix<cub->amx.n;ix++){ \
		    {a} }} /* loop in x-domain */
#define KOOP(a) for(iy=0;iy<ssr->byy.n;iy++){ \
	        for(ix=0;ix<ssr->bxx.n;ix++){ \
		    {a} }} /* loop in k-domain */

#define INDEX(x,a) 0.5+(x-a.o)/a.d;
#define BOUND(i,n) (i<0) ? 0 : ( (i>n-1) ? n-1 : i );
#define  KMAP(i,n) (i<n/2.) ? (i+n/2.) : (i-n/2.);
#define X2K(a,b,p) b.n=a.n+p; b.d=2.0*SF_PI/(b.n*a.d); b.o=(1==b.n)?0:-SF_PI/a.d;

/*------------------------------------------------------------*/
ssr3d ssr3_init( cub3d cub,
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

    /*------------------------------------------------------------*/
    ssr3d ssr;
    ssr = (ssr3d) sf_alloc(1,sizeof(*ssr));

    X2K(cub->amx,ssr->bxx,px);
    X2K(cub->amy,ssr->byy,py);

    /* allocate K-domain storage */
    ssr->wk = sf_complexalloc3 (ssr->bxx.n,ssr->byy.n,cub->ompnth);
    ssr->pk = sf_complexalloc3 (ssr->bxx.n,ssr->byy.n,cub->ompnth);

    /* allocate X-domain storage */
    ssr->wt = sf_floatalloc3   (cub->amx.n,cub->amy.n,cub->ompnth);

    ssr->kk = sf_floatalloc2(ssr->bxx.n,ssr->byy.n);
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
    ssr->lx = sf_intalloc(cub->amx.n);
    ssr->ly = sf_intalloc(cub->amy.n);

    for (iy=0; iy<cub->amy.n; iy++) {
	yy = cub->amy.o + iy*cub->amy.d;
	ily    = INDEX( yy,cub->aly);
	ssr->ly[iy] = BOUND(ily,cub->aly.n); /* x-line index */
    }
    for (ix=0; ix<cub->amx.n; ix++) {
	xx = cub->amx.o + ix*cub->amx.d;
	ilx    = INDEX( xx,cub->alx);
	ssr->lx[ix] = BOUND(ilx,cub->alx.n); /* i-line index */
    }

    /* initialize FFT */
    ssr->f2d = ompfft2_init(cub,ssr->bxx.n,ssr->byy.n);

    ssr->dsmax2 = dsmax*dsmax;
    ssr->dsmax2*= ssr->dsmax2;

    return ssr;
}

/*------------------------------------------------------------*/
void ssr3_close(ssr3d ssr)
/*< free allocated storage >*/
{
    ompfft2_close(ssr->f2d);

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
    cub3d          cub,
    ssr3d          ssr,
    tap3d          tap,
    slo3d          slo,
    int            imz,
    int         ompith
    )
/*< Wavefield extrapolation by SSF >*/
{
    sf_complex w2,co,cc;
    int ix,iy,jr;
    float s,d;
    
#ifdef SF_HAS_COMPLEX_H
    w2 = w*w;
    /* w-x part 1 */
    LOOP( s = 0.5 * slo->so[ompith][ ssr->ly[iy] ][ ssr->lx[ix] ];
	  wx[iy][ix] *= cexpf(-w*s*cub->amz.d); );
#else
    w2 = sf_cmul(w,w);
    /* w-x part 1 */
    LOOP( s = 0.5 * slo->so[ompith][ ssr->ly[iy] ][ ssr->lx[ix] ];
	  wx[iy][ix] = sf_cmul( wx[iy][ix],cexpf(sf_crmul(w,-s*cub->amz.d))); );
#endif

    /* FFT */
    KOOP( ssr->pk[ompith][iy][ix] = sf_cmplx(0.,0.); );
    LOOP( ssr->pk[ompith][iy][ix] = wx[iy][ix]; );
    ompfft2(false,(kiss_fft_cpx**) ssr->pk[ompith],ompith,ssr->f2d);

    LOOP( wx[iy][ix] = sf_cmplx(0.,0.);
	  ssr->wt[ompith][iy][ix] = 0; );
    for (jr=0; jr<slo->nr[imz]; jr++) {
	/* w-k */
#ifdef SF_HAS_COMPLEX_H
	co =       csqrtf(w2 * slo->sm[imz][jr]);
	KOOP( cc = csqrtf(w2 * slo->sm[imz][jr] + ssr->kk[iy][ix]);
	      ssr->wk[ompith][iy][ix] = 
	      ssr->pk[ompith][iy][ix] * cexpf((co-cc)*cub->amz.d); 
	    );
#else
	co =       csqrtf(sf_crmul(w2,slo->sm[imz][jr]));
	KOOP( cc = csqrtf(sf_cadd(sf_crmul(w2,slo->sm[imz][jr]),
				  sf_cmplx(ssr->kk[iy][ix],0.)));
	      ssr->wk[ompith][iy][ix] = 
	      sf_cmul(ssr->pk[ompith][iy][ix],
		      cexpf(sf_crmul(sf_csub(co,cc),cub->amz.d))); 
	    );
#endif
	
	/* IFT */
	ompfft2(true,(kiss_fft_cpx**) ssr->wk[ompith],ompith,ssr->f2d);
	
	/* accumulate wavefield */
#ifdef SF_HAS_COMPLEX_H
	LOOP( d = fabsf(slo->so[ompith][ ssr->ly[iy] ][ ssr->lx[ix] ] * 
			slo->so[ompith][ ssr->ly[iy] ][ ssr->lx[ix] ] - 
			slo->sm[imz][jr]);
	      d = ssr->dsmax2/(d*d+ssr->dsmax2);
	      wx[iy][ix] += ssr->wk[ompith][iy][ix]*d;
	      ssr->wt[ompith][iy][ix] += d; 
	    );
#else
	LOOP( d = fabsf(slo->so[ompith][ ssr->ly[iy] ][ ssr->lx[ix] ] * 
			slo->so[ompith][ ssr->ly[iy] ][ ssr->lx[ix] ] - 
			slo->sm[imz][jr]);
	      d = ssr->dsmax2/(d*d+ssr->dsmax2);
	      wx[iy][ix] = sf_cadd(wx[iy][ix],
				   sf_crmul(ssr->wk[ompith][iy][ix],d));
	      ssr->wt[ompith][iy][ix] += d; 
	    );
#endif
    }
#ifdef SF_HAS_COMPLEX_H
    LOOP( wx[iy][ix] /= ssr->wt[ompith][iy][ix]; );
    /* w-x part 2 */
    LOOP( s = 0.5 * slo->ss[ompith][ ssr->ly[iy] ][ ssr->lx[ix] ];
	  wx[iy][ix] *= cexpf(-w*s*cub->amz.d); );
#else
    LOOP( wx[iy][ix] = sf_crmul(wx[iy][ix],1./ssr->wt[ompith][iy][ix]); );
    /* w-x part 2 */
    LOOP( s = 0.5 * slo->ss[ompith][ ssr->ly[iy] ][ ssr->lx[ix] ];
	  wx[iy][ix] = sf_cmul(wx[iy][ix],cexpf(sf_crmul(w,-s*cub->amz.d))); );
#endif

    taper2d(wx,tap);
}

/*------------------------------------------------------------*/
void ssr3_sso(
    sf_complex    w /* frequency */,
    sf_complex **wx /* wavefield */,
    cub3d          cub,
    ssr3d          ssr,
    tap3d          tap,
    slo3d          slo,
    int            imz,
    int         ompith
    )
/*< Wavefield extrapolation by SSF w/ one reference slowness >*/
{
    sf_complex w2,co,cc;
    int ix,iy,jr;
    float s;

#ifdef SF_HAS_COMPLEX_H
    w2 = w*w;
    /* w-x part 1 */
    LOOP( s = 0.5 * slo->so[ompith][ ssr->ly[iy] ][ ssr->lx[ix] ];
	  wx[iy][ix] *= cexpf(-w*s*cub->amz.d); );
#else
    w2 = sf_cmul(w,w);
    /* w-x part 1 */
    LOOP( s = 0.5 * slo->so[ompith][ ssr->ly[iy] ][ ssr->lx[ix] ];
	  wx[iy][ix] = sf_cmul(wx[iy][ix],cexpf(sf_crmul(w,-s*cub->amz.d))); );
#endif

    /* FFT */
    KOOP( ssr->pk[ompith][iy][ix] = sf_cmplx(0.,0.); );
    LOOP( ssr->pk[ompith][iy][ix] = wx[iy][ix]; );
    ompfft2(false,(kiss_fft_cpx**) ssr->pk[ompith],ompith,ssr->f2d);

    jr=0;
#ifdef SF_HAS_COMPLEX_H
    co =       csqrtf(w2 * slo->sm[imz][jr]);
    KOOP( cc = csqrtf(w2 * slo->sm[imz][jr] + ssr->kk[iy][ix]);
	  ssr->wk[ompith][iy][ix] = 
	  ssr->pk[ompith][iy][ix] * cexpf((co-cc)*cub->amz.d); 
	);

    KOOP( ssr->wk[ompith][iy][ix] = 
	  ssr->pk[ompith][iy][ix] * cexpf((co)*cub->amz.d); );
#else
    co =       csqrtf(sf_crmul(w2,slo->sm[imz][jr]));
    KOOP( cc = csqrtf(sf_cadd(sf_crmul(w2,slo->sm[imz][jr]),
			      sf_cmplx(ssr->kk[iy][ix],0.)));
	  ssr->wk[ompith][iy][ix] = sf_cmul(ssr->pk[ompith][iy][ix],
			       cexpf(sf_crmul(sf_csub(co,cc),cub->amz.d))); 
	);

    KOOP( ssr->wk[ompith][iy][ix] = sf_cmul(ssr->pk[ompith][iy][ix],cexpf(sf_crmul(co,cub->amz.d))); );
#endif

    /* IFT */
    ompfft2(true,(kiss_fft_cpx**) ssr->wk[ompith],ompith,ssr->f2d);
    LOOP( wx[iy][ix] = ssr->wk[ompith][iy][ix]; );
    
    /* w-x part 2 */
#ifdef SF_HAS_COMPLEX_H
    LOOP( s = 0.5 * slo->ss[ompith][ ssr->ly[iy] ][ ssr->lx[ix] ];
	  wx[iy][ix] *= cexpf(-w*s*cub->amz.d); );
#else
    LOOP( s = 0.5 * slo->ss[ompith][ ssr->ly[iy] ][ ssr->lx[ix] ];
	  wx[iy][ix] = sf_cmul(wx[iy][ix],cexpf(sf_crmul(w,-s*cub->amz.d))); );
#endif

    taper2d(wx,tap);
}

/*------------------------------------------------------------*/
void ssr3_phs(
    sf_complex    w /* frequency */,
    sf_complex **wx /* wavefield */,
    cub3d          cub,
    ssr3d          ssr,
    tap3d          tap,
    slo3d          slo,
    int            imz,
    int         ompith
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
    ompfft2(false,(kiss_fft_cpx**) ssr->pk[ompith],ompith,ssr->f2d);

    LOOP( wx[iy][ix] = sf_cmplx(0.,0.);
	  ssr->wt[ompith][iy][ix] = 0; );
    for (jr=0; jr<slo->nr[imz]; jr++) {
	/* w-k */
#ifdef SF_HAS_COMPLEX_H
	KOOP( cc = csqrtf(w2 * slo->sm[imz][jr] + ssr->kk[iy][ix]);
	      ssr->wk[ompith][iy][ix] = 
	      ssr->pk[ompith][iy][ix] * cexpf((-cc)*cub->amz.d); 
	    );
#else
	KOOP( cc = csqrtf(sf_cadd(sf_crmul(w2,slo->sm[imz][jr]),
				  sf_cmplx(ssr->kk[iy][ix],0.)));
	      ssr->wk[ompith][iy][ix] = 
	      sf_cmul(ssr->pk[ompith][iy][ix],cexpf(sf_crmul(cc,-cub->amz.d))); 
	    );
#endif

	/* IFT */
	ompfft2(true,(kiss_fft_cpx**) ssr->wk[ompith],ompith,ssr->f2d);

	/* accumulate wavefield */
#ifdef SF_HAS_COMPLEX_H
	LOOP( d = fabsf(slo->so[ompith][ ssr->ly[iy] ][ ssr->lx[ix] ] * 
			slo->so[ompith][ ssr->ly[iy] ][ ssr->lx[ix] ] - slo->sm[imz][jr]);
	      d = ssr->dsmax2/(d*d+ssr->dsmax2);
	      wx[iy][ix] += ssr->wk[ompith][iy][ix]*d;
	      ssr->wt[ompith][iy][ix] += d; );
#else
	LOOP( d = fabsf(slo->so[ompith][ ssr->ly[iy] ][ ssr->lx[ix] ] * 
			slo->so[ompith][ ssr->ly[iy] ][ ssr->lx[ix] ] - slo->sm[imz][jr]);
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
    cub3d          cub,
    ssr3d          ssr,
    tap3d          tap,
    slo3d          slo,
    int            imz,
    int         ompith
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
    ompfft2(false,(kiss_fft_cpx**) ssr->pk[ompith],ompith,ssr->f2d);

    jr=0;
#ifdef SF_HAS_COMPLEX_H
    KOOP(
	cc = csqrtf(w2 * slo->sm[imz][jr] + ssr->kk[iy][ix]);
	
	ssr->wk[ompith][iy][ix] = 
	ssr->pk[ompith][iy][ix] * cexpf(-cc*cub->amz.d); 
	);
#else
    KOOP(
	cc = csqrtf(sf_cadd(sf_crmul(w2,slo->sm[imz][jr]),
			    sf_cmplx(ssr->kk[iy][ix],0.)));
	
	ssr->wk[ompith][iy][ix] = 
	sf_cmul(ssr->pk[ompith][iy][ix],cexpf(sf_crmul(cc,-cub->amz.d))); 
	);
#endif
    
    /* IFT */
    ompfft2( true,(kiss_fft_cpx**) ssr->wk[ompith],ompith,ssr->f2d);    
    LOOP( wx[iy][ix] = ssr->wk[ompith][iy][ix]; );
    
    taper2d(wx,tap);
    
}

