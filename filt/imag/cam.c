/* 3-D CAM migration/modeling using extended split-step */
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

#include "fft3.h"
#include "taper.h"

#include "slice.h"
/*^*/

#define LOOP(a) for(ihx=0;ihx<ahx.n;ihx++){ \
                for(imy=0;imy<amy.n;imy++){ \
                for(imx=0;imx<amx.n;imx++){ {a} }}}
#define KOOP(a) for(ihx=0;ihx<bhx.n;ihx++){ \
                for(imy=0;imy<bmy.n;imy++){ \
		for(imx=0;imx<bmx.n;imx++){ {a} }}}
#define SOOP(a) for(ily=0;ily<aly.n;ily++){ \
                for(ilx=0;ilx<alx.n;ilx++){ {a} }}

#define INDEX(x,a) 0.5+(x-a.o)/a.d;
#define BOUND(i,n) (i<0) ? 0 : ( (i>n-1) ? n-1 : i );
#define  KMAP(i,n) (i<n/2.) ? (i+n/2.) : (i-n/2.);
#define X2K(a,b,p) b.n=a.n+p; b.d=2.0*SF_PI/(b.n*a.d); b.o=(1==b.n)?0:-SF_PI/a.d;

static sf_axa az,amx,amy,ahx;
static sf_axa    bmx,bmy,bhx;
static sf_axa    alx,aly;
static float dsmax2;

static float          **ksx;/* source   wavenumber  */
static float          **krx;/* receiver wavenumber  */
static int             *jx; /* i-line index */
static int             *jy; /* x-line index */
static int            **is; /* source   index */
static int            **ir; /* receiver index */
static sf_complex ***pk; /* wavefield k */
static sf_complex ***wk; /* wavefield k */
static float         ***wt; /* interpolation weight */

/*------------------------------------------------------------*/
void cam_init(
    sf_axis az_,
    sf_axis amx_,
    sf_axis amy_,
    sf_axis ahx_,
    sf_axis alx_,
    sf_axis aly_,
    int pmx,
    int pmy,
    int phx,
    int tmx,
    int tmy,
    int thx,
    float dsmax
    )
/*< initialize >*/
{
    int   imy, imx, ihx, ilx, ily;
    int        jmx, jhx;
    float  my,  mx,  hx,    k;
    float      kmx, khx;

    az = sf_nod(az_);
    amx= sf_nod(amx_);
    amy= sf_nod(amy_);
    alx= sf_nod(alx_);
    aly= sf_nod(aly_);
    ahx= sf_nod(ahx_);
    
    /* construct K-domain axes */
    X2K(amx,bmx,pmx);
    X2K(amy,bmy,pmy);
    X2K(ahx,bhx,phx);
    fft3_init(bmx.n,bmy.n,bhx.n);

    /* precompute wavenumbers */
    ksx= sf_floatalloc2(bmx.n,bhx.n);
    krx= sf_floatalloc2(bmx.n,bhx.n);
    for (imx=0; imx<bmx.n; imx++) {
	jmx = KMAP(imx,bmx.n);
	kmx = bmx.o + jmx*bmx.d;
	
	for (ihx=0; ihx<bhx.n; ihx++) {
	    jhx = KMAP(ihx,bhx.n);
	    khx = bhx.o + jhx*bhx.d;
	    
	    k = 0.5*(kmx-khx);
	    ksx[ihx][imx] = k*k; /* ksx^2 */
	    
	    k = 0.5*(kmx+khx);
	    krx[ihx][imx] = k*k; /* krx^2 */
	}
    }
    
    /* precompute indices */
    jx = sf_intalloc(amx.n);
    jy = sf_intalloc(amy.n);
    is = sf_intalloc2(amx.n,ahx.n);  /* source   index */
    ir = sf_intalloc2(amx.n,ahx.n);  /* receiver index */

    for (imy=0; imy<amy.n; imy++) {
	my = amy.o + imy*amy.d;
	ily     = INDEX( my,aly);
	jy[imy] = BOUND(ily,aly.n);            /* x-line index */
    }
    for (imx=0; imx<amx.n; imx++) {
	mx = amx.o + imx*amx.d;
	ilx     = INDEX( mx,alx);
	jx[imx] = BOUND(ilx,alx.n);            /* i-line index */
	
	for (ihx=0; ihx<ahx.n; ihx++) {
	    hx = ahx.o + ihx*ahx.d;
	    
	    ilx          = INDEX(mx-hx,alx);
	    is[ihx][imx] = BOUND(  ilx,alx.n); /* source index */
	    
	    ilx          = INDEX(mx+hx,alx);
	    ir[ihx][imx] = BOUND(  ilx,alx.n); /* receiver index */
	}
    }

    /* precompute taper */
    taper3_init(ahx.n,amy.n,amx.n,
		SF_MIN(thx,ahx.n-1),
		SF_MIN(tmy,amy.n-1),
		SF_MIN(tmx,amx.n-1),
		false,true,true);

    /* allocate K-domain storage */
    wk = sf_complexalloc3 (bmx.n,bmy.n,bhx.n);
    pk = sf_complexalloc3 (bmx.n,bmy.n,bhx.n);

    /* allocate X-domain storage */
    wt = sf_floatalloc3   (amx.n,amy.n,ahx.n);

    dsmax2 = dsmax*dsmax;
    dsmax2*= dsmax2;
}

/*------------------------------------------------------------*/

void cam_close(void)
/*< free allocated storage >*/
{
    fft3_close();
    taper3_close();

    ;           free( *ksx);free(ksx);
    ;           free( *krx);free(krx);
    ;           free( *is); free( is);
    ;           free( *ir); free( ir);
    ;                       free( jx);
    ;                       free( jy);

    free(**pk); free( *pk); free( pk);
    free(**wk); free( *wk); free( wk);
    free(**wt); free( *wt); free( wt);
}

/*------------------------------------------------------------*/

void cam_ssf(
    sf_complex      w /* frequency */,
    sf_complex  ***wx /* wavefield */,
    float         **so /* slowness  */, 
    float         **ss /* slowness  */,
    int             nr /* nr. of ref slo */,
    float          *sm /* ref slo squared */
    )
/*< Wavefield extrapolation by SSF >*/
{
    int imy,imx,ihx,jmy,js,jr;
    float  s, kmy, d, dsc,drc;
    sf_complex co, w2, cs, cr, khy, kss, krr;
    
    co = sf_cmplx(0,0);
#ifdef SF_HAS_COMPLEX_H
    w2 = w*w;
#else
    w2 = sf_cmul(w,w);
#endif
    
#ifdef SF_HAS_COMPLEX_H
    LOOP( s = so[ jy[imy] ][ is[ihx][imx] ] 
	  +   so[ jy[imy] ][ ir[ihx][imx] ];
	  wx[ihx][imy][imx] *= cexpf(-w*s* az.d/2); );
#else
    LOOP( s = so[ jy[imy] ][ is[ihx][imx] ] 
	  +   so[ jy[imy] ][ ir[ihx][imx] ];
	  wx[ihx][imy][imx] = sf_cmul(wx[ihx][imy][imx], 
				      cexpf(sf_crmul(w,-s*az.d/2))); );
#endif

     /* FFT */
    KOOP( pk[ihx][imy][imx] = sf_cmplx(0.,0.); );
    LOOP( pk[ihx][imy][imx] = wx[ihx][imy][imx]; );
    fft3(false,(kiss_fft_cpx***) pk);

    LOOP( wx[ihx][imy][imx] = sf_cmplx(0.,0.);
	  wt[ihx][imy][imx] = 0.; );
    for (js=0; js<nr; js++) {
	for (jr=0; jr<nr; jr++) {
	    
	    /* w-k phase shift */
#ifdef SF_HAS_COMPLEX_H
	    co  = csqrtf(w2*sm[js]) 
		+ csqrtf(w2*sm[jr]);
	    KOOP( jmy = KMAP(imy,bmy.n);
		  kmy = bmy.o + jmy*bmy.d; 
		  cs  = csqrtf(w2*sm[js] + ksx[ihx][imx]);
		  cr  = csqrtf(w2*sm[jr] + krx[ihx][imx]);
		  khy = kmy*(cr-cs)/(cr+cs);
		  kss = 0.5*(kmy-khy);
		  krr = 0.5*(kmy+khy);
		  kss = kss*kss + ksx[ihx][imx];
		  krr = krr*krr + krx[ihx][imx];
		  cs  = csqrtf(w2*sm[js] + kss);
		  cr  = csqrtf(w2*sm[jr] + krr);
		  wk[ihx][imy][imx] = 
		  pk[ihx][imy][imx] * cexpf((co-cs-cr)*az.d);
		);
#else
	    co  = sf_cadd(
		csqrtf(sf_crmul(w2,sm[js])),
		csqrtf(sf_crmul(w2,sm[jr])));
	    KOOP( jmy = KMAP(imy,bmy.n);
		  kmy = bmy.o + jmy*bmy.d; 
		  cs  = csqrtf(sf_cadd(sf_crmul(w2,sm[js]),
				       sf_cmplx(ksx[ihx][imx],0.)));
		  cr  = csqrtf(sf_cadd(sf_crmul(w2,sm[jr]),
				       sf_cmplx(krx[ihx][imx],0.)));
		  khy = sf_crmul(sf_cdiv(sf_csub(cr,cs),
					 sf_cadd(cr,cs)),kmy);
		  kss = sf_crmul(sf_csub(sf_cmplx(kmy,0.),khy),0.5);
		  krr = sf_crmul(sf_cadd(sf_cmplx(kmy,0.),khy),0.5);
		  kss = sf_cadd(sf_cmul(kss,kss),sf_cmplx(ksx[ihx][imx],0.));
		  krr = sf_cadd(sf_cmul(krr,krr),sf_cmplx(krx[ihx][imx],0.));
		  cs  = csqrtf(sf_cadd(sf_crmul(w2,sm[js]),kss));
		  cr  = csqrtf(sf_cadd(sf_crmul(w2,sm[jr]),krr));
		  wk[ihx][imy][imx] = sf_cmul(
		      pk[ihx][imy][imx],
		      cexpf(sf_crmul(sf_csub(co,sf_cadd(cs,cr)),az.d)));
		);
#endif
	    
	    /* IFT */
	    fft3(true,(kiss_fft_cpx***) wk);
	    
	    /* accumulate wavefield */
#ifdef SF_HAS_COMPLEX_H
	    LOOP( 
		dsc = fabsf( so[ jy[imy] ][ is[ihx][imx] ] *
			     so[ jy[imy] ][ is[ihx][imx] ] - sm[js]);
		drc = fabsf( so[ jy[imy] ][ ir[ihx][imx] ] *
			     so[ jy[imy] ][ ir[ihx][imx] ] - sm[jr]);
		d = hypotf(dsc,drc);
		d = dsmax2/(d*d+dsmax2);
		wx[ihx][imy][imx] += wk[ihx][imy][imx]*d;
		wt[ihx][imy][imx] += d; );
#else
	    LOOP( 
		dsc = fabsf( so[ jy[imy] ][ is[ihx][imx] ] *
			     so[ jy[imy] ][ is[ihx][imx] ] - sm[js]);
		drc = fabsf( so[ jy[imy] ][ ir[ihx][imx] ] *
			     so[ jy[imy] ][ ir[ihx][imx] ] - sm[jr]);
		d = hypotf(dsc,drc);
		d = dsmax2/(d*d+dsmax2);
		wx[ihx][imy][imx] = sf_cadd(wx[ihx][imy][imx],
					    sf_crmul(wk[ihx][imy][imx],d));
		wt[ihx][imy][imx] += d; );
#endif
	} /* jr loop */
    } /* js loop */
#ifdef SF_HAS_COMPLEX_H
    LOOP( wx[ihx][imy][imx] /= wt[ihx][imy][imx]; );
    LOOP( s = ss[ jy[imy] ][ is[ihx][imx] ] 
	  +   ss[ jy[imy] ][ ir[ihx][imx] ];
	  wx[ihx][imy][imx] *= cexpf(-w*s* az.d/2); );
#else
    LOOP( wx[ihx][imy][imx] = sf_crmul(wx[ihx][imy][imx],
				       1./wt[ihx][imy][imx]); );
    LOOP( s = ss[ jy[imy] ][ is[ihx][imx] ] 
	  +   ss[ jy[imy] ][ ir[ihx][imx] ];
	  wx[ihx][imy][imx] = sf_cmul(wx[ihx][imy][imx],
				      cexpf(sf_crmul(w,-s* az.d/2))); );
#endif

    taper3(wx);
}
