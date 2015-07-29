/* 3-D CAM using extended SSF */

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

#define LOOP(a) for(ihx=0;ihx<cub->ahx.n;ihx++){ \
                for(imy=0;imy<cub->amy.n;imy++){ \
                for(imx=0;imx<cub->amx.n;imx++){ \
		    {a} }}} /* loop in x-domain */
#define KOOP(a) for(ihx=0;ihx<cam->bhx.n;ihx++){ \
                for(imy=0;imy<cam->bmy.n;imy++){ \
		for(imx=0;imx<cam->bmx.n;imx++){ \
		    {a} }}} /* loop in k-domain */

#define INDEX(x,a) 0.5+(x-a.o)/a.d;
#define BOUND(i,n) (i<0) ? 0 : ( (i>n-1) ? n-1 : i );
#define  KMAP(i,n) (i<n/2.) ? (i+n/2.) : (i-n/2.);
#define X2K(a,b,p) b.n=a.n+p; b.d=2.0*SF_PI/(b.n*a.d); b.o=(1==b.n)?0:-SF_PI/a.d;

/*------------------------------------------------------------*/
cam3d cam3_init(cub3d cub,
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
    int   imy, imx, ihx;
    float      kmx, khx;
    int        jmx, jhx;

    int   ilx, ily;

    float  my,  mx,  hx,    k;

    /*------------------------------------------------------------*/
    cam3d cam;
    cam = (cam3d) sf_alloc(1,sizeof(*cam));

    X2K(cub->amx,cam->bmx,pmx);
    X2K(cub->amy,cam->bmy,pmy);
    X2K(cub->ahx,cam->bhx,phx);

    /* allocate K-domain storage */
    cam->wk = sf_complexalloc4 (cam->bmx.n,cam->bmy.n,cam->bhx.n,cub->ompnth);
    cam->pk = sf_complexalloc4 (cam->bmx.n,cam->bmy.n,cam->bhx.n,cub->ompnth);

    /* allocate X-domain storage */
    cam->wt = sf_floatalloc4   (cub->amx.n,cub->amy.n,cub->ahx.n,cub->ompnth);

    cam->ksx = sf_floatalloc2(cam->bmx.n,cam->bhx.n);
    cam->krx = sf_floatalloc2(cam->bmx.n,cam->bhx.n);
    for (imx=0; imx<cam->bmx.n; imx++) {
	jmx = KMAP(imx,cam->bmx.n);
	kmx = cam->bmx.o + jmx*cam->bmx.d;
	
	for (ihx=0; ihx<cam->bhx.n; ihx++) {
	    jhx = KMAP(ihx,cam->bhx.n);
	    khx = cam->bhx.o + jhx*cam->bhx.d;
	    
	    k = 0.5*(kmx-khx);
	    cam->ksx[ihx][imx] = k*k; /* ksx^2 */
	    
	    k = 0.5*(kmx+khx);
	    cam->krx[ihx][imx] = k*k; /* krx^2 */
	}
    }

    /* precompute indices */
    cam->jx = sf_intalloc(cub->amx.n);
    cam->jy = sf_intalloc(cub->amy.n);
    cam->is = sf_intalloc2(cub->amx.n,cub->ahx.n);  /* source   index */
    cam->ir = sf_intalloc2(cub->amx.n,cub->ahx.n);  /* receiver index */

    for (imy=0; imy<cub->amy.n; imy++) {
	my = cub->amy.o + imy*cub->amy.d;
	ily          = INDEX( my,cub->aly);
	cam->jy[imy] = BOUND(ily,cub->aly.n);            /* x-line index */
    }
    for (imx=0; imx<cub->amx.n; imx++) {
	mx = cub->amx.o + imx*cub->amx.d;
	ilx          = INDEX( mx,cub->alx);
	cam->jx[imx] = BOUND(ilx,cub->alx.n);            /* i-line index */
	
	for (ihx=0; ihx<cub->ahx.n; ihx++) {
	    hx = cub->ahx.o + ihx*cub->ahx.d;
	    
	    ilx               = INDEX(mx-hx,cub->alx);
	    cam->is[ihx][imx] = BOUND(  ilx,cub->alx.n); /* source index */
	    
	    ilx               = INDEX(mx+hx,cub->alx);
	    cam->ir[ihx][imx] = BOUND(  ilx,cub->alx.n); /* receiver index */
	}
    }

    /* initialize FFT */
    cam->f3d = ompfft3_init(cub,cam->bmx.n,cam->bmy.n,cam->bhx.n);
    
    cam->dsmax2 = dsmax*dsmax;
    cam->dsmax2*= cam->dsmax2;

    return cam;
}

/*------------------------------------------------------------*/

void cam3_close(cam3d cam)
/*< free allocated storage >*/
{
    ompfft3_close(cam->f3d);

    ;           free( *cam->ksx);free(cam->ksx);
    ;           free( *cam->krx);free(cam->krx);
    ;           free( *cam->is); free(cam->is);
    ;           free( *cam->ir); free(cam->ir);
    ;                       free( cam->jx);
    ;                       free( cam->jy);

    free(**cam->pk); free( *cam->pk); free( cam->pk);
    free(**cam->wk); free( *cam->wk); free( cam->wk);
    free(**cam->wt); free( *cam->wt); free( cam->wt);
}

/*------------------------------------------------------------*/
void cam3_ssf(
    sf_complex       w /* frequency */,
    sf_complex   ***wx /* wavefield */,
    cub3d          cub,
    cam3d          cam,
    tap3d          tap,
    slo3d          slo,
    int            imz,
    int         ompith
    )
/*< Wavefield extrapolation by SSF >*/
{
    sf_complex co,cs,cr,w2,khy,kss,krr;
    float s, kmy, d,dsc,drc;
    int imy,imx,ihx,jmy, js,jr;

    co = sf_cmplx(0,0);
#ifdef SF_HAS_COMPLEX_H
    w2 = w*w;
#else
    w2 = sf_cmul(w,w);
#endif
    
#ifdef SF_HAS_COMPLEX_H
    LOOP( s = slo->so[ompith][ cam->jy[imy] ][ cam->is[ihx][imx] ] 
	  +   slo->so[ompith][ cam->jy[imy] ][ cam->ir[ihx][imx] ];
	  wx[ihx][imy][imx] *= cexpf(-w*s* cub->amz.d/2); );
#else
    LOOP( s = slo->so[ompith][ cam->jy[imy] ][ cam->is[ihx][imx] ] 
	  +   slo->so[ompith][ cam->jy[imy] ][ cam->ir[ihx][imx] ];
	  wx[ihx][imy][imx] = sf_cmul(wx[ihx][imy][imx], 
				      cexpf(sf_crmul(w,-s*cub->amz.d/2))); );
#endif

    /* FFT */
    KOOP( cam->pk[ompith][ihx][imy][imx] = sf_cmplx(0.,0.); );
    LOOP( cam->pk[ompith][ihx][imy][imx] = wx[ihx][imy][imx]; );
    ompfft3(false,(kiss_fft_cpx***) cam->pk[ompith],ompith,cam->f3d);

    LOOP( wx[ihx][imy][imx] = sf_cmplx(0.,0.);
	  cam->wt[ompith][ihx][imy][imx] = 0.; );

    for (js=0; js<slo->nr[imz]; js++) {
	for (jr=0; jr<slo->nr[imz]; jr++) {
	    
	    /* w-k phase shift */
#ifdef SF_HAS_COMPLEX_H
	    co  = csqrtf(w2*slo->sm[imz][js]) 
		+ csqrtf(w2*slo->sm[imz][jr]);

	    KOOP( jmy = KMAP(imy,cam->bmy.n);
		  kmy = cam->bmy.o + jmy*cam->bmy.d; 
		  cs  = csqrtf(w2*slo->sm[imz][js] + cam->ksx[ihx][imx]);
		  cr  = csqrtf(w2*slo->sm[imz][jr] + cam->krx[ihx][imx]);
		  khy = kmy*(cr-cs)/(cr+cs);
		  kss = 0.5*(kmy-khy);
		  krr = 0.5*(kmy+khy);
		  kss = kss*kss + cam->ksx[ihx][imx];
		  krr = krr*krr + cam->krx[ihx][imx];
		  cs  = csqrtf(w2*slo->sm[imz][js] + kss);
		  cr  = csqrtf(w2*slo->sm[imz][jr] + krr);
		  cam->wk[ompith][ihx][imy][imx] = 
		  cam->pk[ompith][ihx][imy][imx] * cexpf((co-cs-cr)*cub->amz.d);
		);
#else
	    co  = sf_cadd(
		csqrtf(sf_crmul(w2,slo->sm[imz][js])),
		csqrtf(sf_crmul(w2,slo->sm[imz][jr])));
	    KOOP( jmy = KMAP(imy,cam->bmy.n);
		  kmy = cam->bmy.o + jmy*cam->bmy.d; 
		  cs  = csqrtf(sf_cadd(sf_crmul(w2,slo->sm[imz][js]),
				       sf_cmplx(cam->ksx[ihx][imx],0.)));
		  cr  = csqrtf(sf_cadd(sf_crmul(w2,slo->sm[imz][jr]),
				       sf_cmplx(cam->krx[ihx][imx],0.)));
		  khy = sf_crmul(sf_cdiv(sf_csub(cr,cs),
					 sf_cadd(cr,cs)),kmy);
		  kss = sf_crmul(sf_csub(sf_cmplx(kmy,0.),khy),0.5);
		  krr = sf_crmul(sf_cadd(sf_cmplx(kmy,0.),khy),0.5);
		  kss = sf_cadd(sf_cmul(kss,kss),sf_cmplx(cam->ksx[ihx][imx],0.));
		  krr = sf_cadd(sf_cmul(krr,krr),sf_cmplx(cam->krx[ihx][imx],0.));
		  cs  = csqrtf(sf_cadd(sf_crmul(w2,slo->sm[imz][js]),kss));
		  cr  = csqrtf(sf_cadd(sf_crmul(w2,slo->sm[imz][jr]),krr));
		  cam->wk[ompith][ihx][imy][imx] = sf_cmul(
		  cam->pk[ompith][ihx][imy][imx],
		  cexpf(sf_crmul(sf_csub(co,sf_cadd(cs,cr)),cub->amz.d)));
		);
#endif
	    
	    /* IFT */
	    ompfft3(true,(kiss_fft_cpx***) cam->wk[ompith],ompith,cam->f3d);

	    /* accumulate wavefield */
#ifdef SF_HAS_COMPLEX_H
	    LOOP( 
		dsc = fabsf( slo->so[ompith][ cam->jy[imy] ][ cam->is[ihx][imx] ] *
			     slo->so[ompith][ cam->jy[imy] ][ cam->is[ihx][imx] ] - slo->sm[imz][js]);
		drc = fabsf( slo->so[ompith][ cam->jy[imy] ][ cam->ir[ihx][imx] ] *
			     slo->so[ompith][ cam->jy[imy] ][ cam->ir[ihx][imx] ] - slo->sm[imz][jr]);
		d = hypotf(dsc,drc);
		d = cam->dsmax2/(d*d+cam->dsmax2);
		wx[ihx][imy][imx] += cam->wk[ompith][ihx][imy][imx]*d;
		cam->wt[ompith][ihx][imy][imx] += d; );
#else
	    LOOP( 
		dsc = fabsf( slo->so[ompith][ cam->jy[imy] ][ cam->is[ihx][imx] ] *
			     slo->so[ompith][ cam->jy[imy] ][ cam->is[ihx][imx] ] - slo->sm[imz][js]);
		drc = fabsf( slo->so[ompith][ cam->jy[imy] ][ cam->ir[ihx][imx] ] *
			     slo->so[ompith][ cam->jy[imy] ][ cam->ir[ihx][imx] ] - slo->sm[imz][jr]);
		d = hypotf(dsc,drc);
		d = cam->dsmax2/(d*d+cam->dsmax2);
		wx[ihx][imy][imx] = sf_cadd(wx[ihx][imy][imx],
					    sf_crmul(cam->wk[ompith][ihx][imy][imx],d));
		cam->wt[ompith][ihx][imy][imx] += d; );
#endif

	} /* jr loop */
    } /* js loop */  

#ifdef SF_HAS_COMPLEX_H
    LOOP( wx[ihx][imy][imx] /= cam->wt[ompith][ihx][imy][imx]; );
    LOOP( s = slo->ss[ompith][ cam->jy[imy] ][ cam->is[ihx][imx] ] 
	  +   slo->ss[ompith][ cam->jy[imy] ][ cam->ir[ihx][imx] ];
	  wx[ihx][imy][imx] *= cexpf(-w*s* cub->amz.d/2); );
#else
    LOOP( wx[ihx][imy][imx] = sf_crmul(wx[ihx][imy][imx],
				       1./cam->wt[ompith][ihx][imy][imx]); );
    LOOP( s = slo->ss[ompith][ cam->jy[imy] ][ cam->is[ihx][imx] ] 
	  +   slo->ss[ompith][ cam->jy[imy] ][ cam->ir[ihx][imx] ];
	  wx[ihx][imy][imx] = sf_cmul(wx[ihx][imy][imx],
				      cexpf(sf_crmul(w,-s* cub->amz.d/2))); );
#endif
    
    taper3d(wx,tap);

}
