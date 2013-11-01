/* 3-D SSR using extended SSF */

/*
  Copyright (C) 2010 Colorado School of Mines
  
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

#include "weiutl.h"
/*^*/

#include "weifft.h"

/*------------------------------------------------------------*/

#define LOOP(a) for(iy=0;iy<sf_n(cub->amy);iy++){ \
	        for(ix=0;ix<sf_n(cub->amx);ix++){ \
		    {a} }} /* loop in x-domain */
#define KOOP(a) for(iy=0;iy<sf_n(ssr->byy);iy++){ \
  	        for(ix=0;ix<sf_n(ssr->bxx);ix++){ \
		    {a} }} /* loop in k-domain */

#define INDEX(x,a) 0.5+(x-sf_o(a))/sf_d(a);
#define BOUND(i,n) (i<0) ? 0 : ( (i>n-1) ? n-1 : i );
#define  KMAP(i,n) (i<n/2.) ? (i+n/2.) : (i-n/2.);
#define X2K(a,b,p) sf_setn(b,sf_n(a)+p);		   \
                   sf_setd(b,2.0*SF_PI/(sf_n(b)*sf_d(a)));		   \
                   sf_seto(b,(1==sf_n(b))?0:-SF_PI/sf_d(a));

/*------------------------------------------------------------*/
weissr3d weissr_init(weicub3d cub,
                     weislo3d slo)
/*< initialize >*/
{
    int   ix, iy, iz, jr;
    float kx, ky;
    int   jx, jy;
    float xx, yy;
    int  ilx,ily;
    int ompith;
    float **sum, d;
    weissr3d ssr;

    d = 0.0;
    sum = sf_floatalloc2(sf_n(cub->amx),sf_n(cub->amy));

    /*------------------------------------------------------------*/
    ssr = (weissr3d) sf_alloc(1,sizeof(*ssr));

    ssr->bxx=sf_maxa(1,0,1);
    ssr->byy=sf_maxa(1,0,1);

    X2K(cub->amx,ssr->bxx,cub->pmx);
    X2K(cub->amy,ssr->byy,cub->pmy);

    /* allocate K-domain storage */
    ssr->wk = sf_complexalloc3 (sf_n(ssr->bxx),sf_n(ssr->byy),cub->ompnth);
    ssr->pk = sf_complexalloc3 (sf_n(ssr->bxx),sf_n(ssr->byy),cub->ompnth);

    /* allocate X-domain storage */
    ssr->wt  = sf_floatalloc3(sf_n(cub->amx),sf_n(cub->amy),cub->ompnth);
    ssr->wt1 = sf_floatalloc4(sf_n(cub->amx),sf_n(cub->amy),sf_n(cub->az),slo->nrmax);

    ssr->kk = sf_floatalloc2(sf_n(ssr->bxx),sf_n(ssr->byy));
    for (iy=0; iy<sf_n(ssr->byy); iy++) {
	jy = KMAP(iy,sf_n(ssr->byy));
	ky = sf_o(ssr->byy) + jy*sf_d(ssr->byy);
	for (ix=0; ix<sf_n(ssr->bxx); ix++) {
	    jx = KMAP(ix,sf_n(ssr->bxx));
	    kx = sf_o(ssr->bxx) + jx*sf_d(ssr->bxx);  
	    ssr->kk[iy][ix] = kx*kx+ky*ky;
	}
    }
    
    /* precompute indices */
    ssr->lx = sf_intalloc(sf_n(cub->amx));
    ssr->ly = sf_intalloc(sf_n(cub->amy));

    /* x-line index */
    for (iy=0; iy<sf_n(cub->amy); iy++) {
	yy = sf_o(cub->amy) + iy*sf_d(cub->amy);
	ily    = INDEX( yy,cub->aly);
	ssr->ly[iy] = BOUND(ily,sf_n(cub->aly));
    }
    /* i-line index */
    for (ix=0; ix<sf_n(cub->amx); ix++) {
	xx = sf_o(cub->amx) + ix*sf_d(cub->amx);
	ilx    = INDEX( xx,cub->alx);
	ssr->lx[ix] = BOUND(ilx,sf_n(cub->alx));
    }

    /* initialize FFT */
    ssr->f2d = (weifft2d*) sf_alloc(cub->ompnth,sizeof(weifft2d));
    for(ompith=0;ompith<cub->ompnth;ompith++) {
	ssr->f2d[ompith] = weifft_init(cub,sf_n(ssr->bxx),sf_n(ssr->byy));
    }

    ssr->dsmax2 = cub->dsmax*cub->dsmax;
    ssr->dsmax2*= ssr->dsmax2;

    for(iz=0; iz<sf_n(cub->az); iz++){
        LOOP( sum[iy][ix] = 0.0; );

        for(jr=0; jr<slo->nr[iz]; jr++){
            LOOP( d = fabsf(slo->s[iz][ ssr->ly[iy] ][ ssr->lx[ix] ] *
                            slo->s[iz][ ssr->ly[iy] ][ ssr->lx[ix] ] -
                            slo->sm[iz][jr]);
                  d = ssr->dsmax2/(d*d+ssr->dsmax2);
                  sum[iy][ix] += d;
                );
        }

        for(jr=0; jr<slo->nr[iz]; jr++){
            LOOP( d = fabsf(slo->s[iz][ ssr->ly[iy] ][ ssr->lx[ix] ] *
                            slo->s[iz][ ssr->ly[iy] ][ ssr->lx[ix] ] -
                            slo->sm[iz][jr]);
                  d = ssr->dsmax2/(d*d+ssr->dsmax2);
                  ssr->wt1[jr][iz][iy][ix] = d/sum[iy][ix];
                )
        }
    }

    free( *sum); free(  sum);

    return ssr;
}

/*------------------------------------------------------------*/
void weissr_close(weicub3d cub,
		  weissr3d ssr)
/*< free allocated storage >*/
{
    int ompith;

    /* close FFT */
    for(ompith=0;ompith<cub->ompnth;ompith++) {
	weifft_close(ssr->f2d[ompith]);
    }
    free(ssr->f2d);

    free( *ssr->kk); free( ssr->kk);
    ;                free( ssr->lx);
    ;                free( ssr->ly);

    free(**ssr->pk); free( *ssr->pk); free( ssr->pk);
    free(**ssr->wk); free( *ssr->wk); free( ssr->wk);
    free(**ssr->wt); free( *ssr->wt); free( ssr->wt);
}

/*------------------------------------------------------------*/
void weissr( sf_complex          w, 
	     sf_complex       **wx,
	     weicub3d          cub,
	     weissr3d          ssr,
	     weislo3d          slo,
	     int                iz,
	     int            ompith
    )
/*< Wavefield extrapolation by SSF >*/
{
    sf_complex w2,co,cc;
    int ix,iy,jr;
    float s,d;
     
    /* w-x part 1 */
#ifdef SF_HAS_COMPLEX_H
    w2 = w*w;
    LOOP( s = 0.5 * slo->s[iz][ ssr->ly[iy] ][ ssr->lx[ix] ];
	  wx[iy][ix] *= cexpf(-w*s*sf_d(cub->az)); );
#else
    w2 = sf_cmul(w,w);
    LOOP( s = 0.5 * slo->s[iz][ ssr->ly[iy] ][ ssr->lx[ix] ];
	  wx[iy][ix] = sf_cmul( wx[iy][ix],cexpf(sf_crmul(w,-s*sf_d(cub->az)))); );
#endif
   /* FFT */
    KOOP( ssr->pk[ompith][iy][ix] = sf_cmplx(0.,0.); );
    LOOP( ssr->pk[ompith][iy][ix] = wx[iy][ix]; );
    weifft(false,
	   (kiss_fft_cpx**) ssr->pk[ompith],
	   ssr->f2d[ompith]);

    /* w-k */
    LOOP( wx[iy][ix] = sf_cmplx(0.,0.);
	  ssr->wt[ompith][iy][ix] = 0; );
    for (jr=0; jr<slo->nr[iz]; jr++) {
#ifdef SF_HAS_COMPLEX_H
	co =       csqrtf(w2 * slo->sm[iz][jr]);
	KOOP( cc = csqrtf(w2 * slo->sm[iz][jr] + ssr->kk[iy][ix]);
	      ssr->wk[ompith][iy][ix] = 
	      ssr->pk[ompith][iy][ix] * cexpf((co-cc)*sf_d(cub->az)); 
	    );
#else
	co =       csqrtf(sf_crmul(w2,slo->sm[iz][jr]));
	KOOP( cc = csqrtf(sf_cadd(sf_crmul(w2,slo->sm[iz][jr]),
				  sf_cmplx(ssr->kk[iy][ix],0.)));
	      ssr->wk[ompith][iy][ix] = 
	      sf_cmul(ssr->pk[ompith][iy][ix],
		      cexpf(sf_crmul(sf_csub(co,cc),sf_d(cub->az)))); 
	    );
#endif
	
	/* IFT */
	weifft(true,
	       (kiss_fft_cpx**) ssr->wk[ompith],
	       ssr->f2d[ompith]);

	/* accumulate wavefield */
#ifdef SF_HAS_COMPLEX_H
	LOOP( d = fabsf(slo->s[iz][ ssr->ly[iy] ][ ssr->lx[ix] ] * 
			slo->s[iz][ ssr->ly[iy] ][ ssr->lx[ix] ] - 
			slo->sm[iz][jr]);
	      d = ssr->dsmax2/(d*d+ssr->dsmax2);
	      wx[iy][ix] += ssr->wk[ompith][iy][ix]*d;
	      ssr->wt[ompith][iy][ix] += d; 
	    );
#else
	LOOP( d = fabsf(slo->s[iz][ ssr->ly[iy] ][ ssr->lx[ix] ] * 
			slo->s[iz][ ssr->ly[iy] ][ ssr->lx[ix] ] - 
			slo->sm[iz][jr]);
	      d = ssr->dsmax2/(d*d+ssr->dsmax2);
	      wx[iy][ix] = sf_cadd(wx[iy][ix],
			   sf_crmul(ssr->wk[ompith][iy][ix],d));
	      ssr->wt[ompith][iy][ix] += d; 
	    );
#endif

    }

    /* w-x part 2 */
#ifdef SF_HAS_COMPLEX_H
    LOOP( wx[iy][ix] /= ssr->wt[ompith][iy][ix]; );
    LOOP( s = 0.5 * slo->s[iz+1][ ssr->ly[iy] ][ ssr->lx[ix] ];
	  wx[iy][ix] *= cexpf(-w*s*sf_d(cub->az)); );
#else
    LOOP( wx[iy][ix] = sf_crmul(wx[iy][ix],1./ssr->wt[ompith][iy][ix]); );
    LOOP( s = 0.5 * slo->s[iz+1][ ssr->ly[iy] ][ ssr->lx[ix] ];
	  wx[iy][ix] = sf_cmul(wx[iy][ix],cexpf(sf_crmul(w,-s*sf_d(cub->az)))); );
#endif

}

/*------------------------------------------------------------*/
void weissr1(sf_complex          w,
             sf_complex       **wx,
             weicub3d          cub,
             weissr3d          ssr,
             weislo3d          slo,
             int                iz,
             int            ompith,
             bool               adj
    )
/*< Wavefield extrapolation by SSF >*/
{
    sf_complex w2,co,cc;
    int ix,iy,jr,sign,flg;
    float s;
    sf_complex **tmp, **pk, **wk;

    sign = adj?+1:-1; /* adj=true:  downward continuation */
    flg = adj?+0:-1; /* adj=true:  downward continuation */

    tmp = sf_complexalloc2(sf_n(cub->amx),sf_n(cub->amy));
    pk  = sf_complexalloc2(sf_n(ssr->bxx),sf_n(ssr->byy));
    wk  = sf_complexalloc2(sf_n(ssr->bxx),sf_n(ssr->byy));

    /* w-x part 1 */
#ifdef SF_HAS_COMPLEX_H
    w2 = w*w;
    LOOP( s = 0.5 * slo->s[iz][ ssr->ly[iy] ][ ssr->lx[ix] ];
          wx[iy][ix] *= cexpf(-w*s*sf_d(cub->az)); );
#else
    w2 = sf_cmul(w,w);
    LOOP( s = 0.5 * slo->s[iz][ ssr->ly[iy] ][ ssr->lx[ix] ];
          wx[iy][ix] = sf_cmul( wx[iy][ix],cexpf(sf_crmul(w,-s*sf_d(cub->az)))); );
#endif

    LOOP( tmp[iy][ix] = wx[iy][ix];
          wx[iy][ix] = sf_cmplx(0.0,0.0); );

    for (jr=0; jr<slo->nr[iz+flg]; jr++) {
#ifdef SF_HAS_COMPLEX_H
        LOOP(
              if(adj) pk[iy][ix] = tmp[iy][ix];
              else    pk[iy][ix] = tmp[iy][ix]*ssr->wt1[jr][iz+flg][iy][ix];
            );
#else
        LOOP(
              if(adj) pk[iy][ix] = tmp[iy][ix];
              else    pk[iy][ix] = sf_crmul(tmp[iy][ix],ssr->wt1[jr][iz+flg][iy][ix]);
            );
#endif
        /* FFT */
        weifft(false, (kiss_fft_cpx**) pk, ssr->f2d[ompith]);

        /* w-k */
#ifdef SF_HAS_COMPLEX_H
        co =       csqrtf(w2 * slo->sm[iz+flg][jr]);
        KOOP( cc = csqrtf(w2 * slo->sm[iz+flg][jr] + ssr->kk[iy][ix]);
              wk[iy][ix] = pk[iy][ix] * cexpf((co-cc)*sf_d(cub->az));
            );
#else
        co =       csqrtf(sf_crmul(w2,slo->sm[iz+flg][jr]));
        KOOP( cc = csqrtf(sf_cadd(sf_crmul(w2,slo->sm[iz+flg][jr]),
                                  sf_cmplx(ssr->kk[iy][ix],0.)));
              wk[iy][ix] = sf_cmul(pk[iy][ix],
				   cexpf(sf_crmul(sf_csub(co,cc),sf_d(cub->az))));
            );
#endif

        /* IFT */
        weifft(true,(kiss_fft_cpx**) wk,ssr->f2d[ompith]);

        /* accumulate wavefield  */
#ifdef SF_HAS_COMPLEX_H
        LOOP(
              if(adj) wx[iy][ix] += wk[iy][ix]*ssr->wt1[jr][iz+flg][iy][ix];
              else    wx[iy][ix] += wk[iy][ix];
            );
#else
        LOOP(
              if(adj) wx[iy][ix] = sf_cadd(wx[iy][ix],sf_crmul(wk[iy][ix],ssr->wt1[jr][iz+flg][iy][ix]));
              else    wx[iy][ix] = sf_cadd(wx[iy][ix],wk[iy][ix]);
            );
#endif
    }

    /* w-x part 2 */
#ifdef SF_HAS_COMPLEX_H
    LOOP( s = 0.5 * slo->s[iz+sign][ ssr->ly[iy] ][ ssr->lx[ix] ];
          wx[iy][ix] *= cexpf(-w*s*sf_d(cub->az)); );
#else
    LOOP( s = 0.5 * slo->s[iz+sign][ ssr->ly[iy] ][ ssr->lx[ix] ];
          wx[iy][ix] = sf_cmul(wx[iy][ix],cexpf(sf_crmul(w,-s*sf_d(cub->az)))); );
#endif

    free( *tmp); free( tmp);
    free(  *wk); free( wk);
    free(  *pk); free( pk);
}
    
