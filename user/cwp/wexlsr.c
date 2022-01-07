/* 3-D linearized SSR */

/*
  Copyright (C) 2009 Colorado School of Mines
  
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

#include "wexfft.h"
#include "wexutl.h"
/*^*/

/*------------------------------------------------------------*/
#define LOOP(a) for(iy=0;iy<cub->amy.n;iy++){ for(ix=0;ix<cub->amx.n;ix++){ {a} }}
#define KOOP(a) for(iy=0;iy<lsr->byy.n;iy++){ for(ix=0;ix<lsr->bxx.n;ix++){ {a} }}

#define KMAP(i,n) (i<n/2.) ? (i+n/2.) : (i-n/2.);
#define X2K(a,b,p) b.n=a.n+p; b.d=2.0*SF_PI/(b.n*a.d); b.o=(1==b.n)?0:-SF_PI/a.d;
#define INDEX(x,a) 0.5+(x-a.o)/a.d;
#define BOUND(i,n) (i<0) ? 0 : ( (i>n-1) ? n-1 : i );

/*------------------------------------------------------------*/

wexlsr3d wexlsr_init(wexcub3d cub,
                      int      px,
                      int      py,
                      float  dsmax
    )
/*< initialize >*/
{
    int   ix, iy;
    float kx, ky;
    int   jx, jy;
    float xx, yy;
    int  ilx,ily;
    int  ompith;

    /*------------------------------------------------------------*/
    wexlsr3d lsr;
    lsr = (wexlsr3d) sf_alloc(1,sizeof(*lsr));

    X2K(cub->amx,lsr->bxx,px);
    X2K(cub->amy,lsr->byy,py);

    /* allocate K-domain storage */
    lsr->wk = sf_complexalloc3(lsr->bxx.n,lsr->byy.n,cub->ompnth);
    lsr->wt = sf_complexalloc3(lsr->bxx.n,lsr->byy.n,cub->ompnth);

    /* precompute wavenumbers */
    lsr->kk = sf_floatalloc2(lsr->bxx.n,lsr->byy.n);
    lsr->kw = sf_floatalloc2(lsr->bxx.n,lsr->byy.n);
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

    /* precompute indices */
    lsr->lx = sf_intalloc(cub->amx.n);
    lsr->ly = sf_intalloc(cub->amy.n);

    for (iy=0; iy<cub->amy.n; iy++) {
        yy = cub->amy.o + iy*cub->amy.d;
        ily    = INDEX( yy,cub->aly);
        lsr->ly[iy] = BOUND(ily,cub->aly.n); /* x-line index */
    }
    for (ix=0; ix<cub->amx.n; ix++) {
        xx = cub->amx.o + ix*cub->amx.d;
        ilx    = INDEX( xx,cub->alx);
        lsr->lx[ix] = BOUND(ilx,cub->alx.n); /* i-line index */
    }

    /* square-root expansion coefficients */
    if(!sf_getint("nsc",&lsr->nsc)) lsr->nsc = 0;
    if(lsr->nsc>5) lsr->nsc=5;

    lsr->csc[0]= 1.;
    lsr->csc[1]= 1./  2.;
    lsr->csc[2]= 3./  8.;
    lsr->csc[3]= 5./ 16.;
    lsr->csc[4]=35./128.;

    /* initialize FFT */
    lsr->f2d = (wexfft2d*) sf_alloc(cub->ompnth,sizeof(wexfft2d));
    for(ompith=0;ompith<cub->ompnth;ompith++) {
        lsr->f2d[ompith] = wexfft_init(cub,lsr->bxx.n,lsr->byy.n);
    }

    lsr->dsmax2 = dsmax*dsmax;
    lsr->dsmax2*= lsr->dsmax2;

    return lsr;

}

/*------------------------------------------------------------*/

void wexlsr_close(wexcub3d cub,
                  wexlsr3d lsr)
/*< free allocated storage >*/
{

    int ompith;

    /* close FFT */
    for(ompith=0;ompith<cub->ompnth;ompith++) {
        wexfft_close(lsr->f2d[ompith]);
    }
    free(lsr->f2d);

    ;           free( lsr->lx);
    ;           free( lsr->ly);

    free(**lsr->wk); free( *lsr->wk); free( lsr->wk);
    free(**lsr->wt); free( *lsr->wt); free( lsr->wt);

                   ; free( *lsr->kk); free( lsr->kk);
                   ; free( *lsr->kw); free( lsr->kw);
}

/*------------------------------------------------------------*/

void kweight(wexcub3d cub,
             wexlsr3d lsr,
             float **bs,
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

    ko  = fabsf(wo) * smin;
    ko *= ko;

    KOOP(
        if( lsr->kk[iy][ix] < ko ) {
            lsr->kw[iy][ix] = 1.;
        } else {
            lsr->kw[iy][ix] = 0.;
        } );
}

/*------------------------------------------------------------*/
void wexlsr_s2w(
    sf_complex    w /* frequency */,
    sf_complex **bw /* background   wavefield */,
    wexcub3d    cub,
    wexlsr3d    lsr,
    wexslo3d    slo,
    sf_complex **pw /* perturbation wavefield */,
    sf_complex **ps /* perturbation slowness  */,
    int          iz,
    int      ompith
    )
/*< linear scattering operator (forward) >*/
{
    int ix,iy,isc;
    float wo;
    sf_complex iwdz;

    wo = cimagf(w);     /* real frequency */
    iwdz = sf_cmplx(0.,-2 * wo * cub->az.d);

    // 0th order term 
#ifdef SF_HAS_COMPLEX_H
    LOOP( pw[iy][ix] = ps[iy][ix] * iwdz * bw[iy][ix]; );
    //LOOP( pw[iy][ix] = ps[iy][ix] * bw[iy][ix]; );
#else
    //LOOP( pw[iy][ix] = sf_cmul(ps[iy][ix],bw[iy][ix]); );
    LOOP( pw[iy][ix] = sf_cmul(ps[iy][ix], sf_cmul(iwdz,bw[iy][ix])); );
#endif

    // higher order terms 
    if(lsr->nsc>0) {
        kweight(cub,lsr,slo->s[iz],wo);     // k-domain weight 

        LOOP( lsr->wt[ompith][iy][ix] = pw[iy][ix]; );

        for( isc=1; isc<=lsr->nsc; isc++) {
            KOOP( lsr->wk[ompith][iy][ix] = sf_cmplx(0.,0.); );
            LOOP( lsr->wk[ompith][iy][ix] = lsr->wt[ompith][iy][ix]; );

            wexfft(false,(kiss_fft_cpx**) lsr->wk[ompith],lsr->f2d[ompith]);


#ifdef SF_HAS_COMPLEX_H
            KOOP( lsr->wk[ompith][iy][ix] *=
                  lsr->kw[iy][ix] * pow(lsr->kk[iy][ix],isc); );
#else
            KOOP( lsr->wk[ompith][iy][ix] = sf_crmul(lsr->wk[ompith][iy][ix],
                                        lsr->kw[iy][ix] * pow(lsr->kk[iy][ix],isc)); );
#endif

            // IFT 
            wexfft(true,(kiss_fft_cpx**) lsr->wk[ompith],lsr->f2d[ompith]);

#ifdef SF_HAS_COMPLEX_H
            LOOP( lsr->wk[ompith][iy][ix] *= pow(wo*slo->s[iz][ lsr->ly[iy] ][ lsr->lx[ix] ],-2*isc); );
            LOOP( pw[iy][ix] +=
                  lsr->wk[ompith][iy][ix] * lsr->csc[isc]; );
#else
             LOOP( lsr->wk[ompith][iy][ix] = sf_crmul(lsr->wk[ompith][iy][ix],
                                         pow(wo*slo->s[iz][ lsr->ly[iy] ][ lsr->lx[ix] ],-2*isc)); );
             LOOP( pw[iy][ix] = sf_cadd(pw[iy][ix],
                                        sf_crmul(lsr->wk[ompith][iy][ix],lsr->csc[isc])); );
#endif
        }
    }  
}

/*------------------------------------------------------------*/
void wexlsr_w2s(
    sf_complex    w /* frequency */,
    sf_complex **bw /* background   wavefield */,
    wexcub3d    cub,
    wexlsr3d    lsr,
    wexslo3d    slo,
    sf_complex **pw /* perturbation wavefield */,
    sf_complex **ps /* perturbation slowness  */,
    int          iz,
    int      ompith
    )
/*< linear scattering operator (adjoint) >*/
{
    int ix,iy,isc;
    float wo;
    sf_complex iwdz;

    wo = cimagf(w);     /* real frequency */
    iwdz = sf_cmplx(0.,-2 * wo * cub->az.d);

    // 0th order term
#ifdef SF_HAS_COMPLEX_H
    //LOOP( ps[iy][ix] = pw[iy][ix] * conjf(bw[iy][ix]); );
    LOOP( ps[iy][ix] = pw[iy][ix] * iwdz * conjf( bw[iy][ix] ); );
#else
    //LOOP( ps[iy][ix] = sf_cmul(pw[iy][ix],conjf(bw[iy][ix])); );
    LOOP( ps[iy][ix] = sf_cmul(pw[iy][ix],sf_cmul(iwdz,conjf(bw[iy][ix])) ); );
#endif

    // higher order terms 
    if(lsr->nsc>0) {
        kweight(cub,lsr,slo->s[iz],wo);     // k-domain weight 

        for( isc=1; isc<=lsr->nsc; isc++) {
            KOOP( lsr->wk[ompith][iy][ix] = sf_cmplx(0.,0.); );
            LOOP( lsr->wk[ompith][iy][ix] = pw[iy][ix]; );
#ifdef SF_HAS_COMPLEX_H
            LOOP( lsr->wk[ompith][iy][ix]*= pow(wo*slo->s[iz][ lsr->ly[iy] ][lsr->lx[ix]],-2*isc); );
#else
            LOOP( lsr->wk[ompith][iy][ix] = sf_crmul(lsr->wk[ompith][iy][ix],
                                        pow(wo*slo->s[iz][lsr->ly[iy]][lsr->lx[ix]],-2*isc)); );
#endif

           // IFT 
           wexfft(true,(kiss_fft_cpx**) lsr->wk[ompith],lsr->f2d[ompith]);


#ifdef SF_HAS_COMPLEX_H
            KOOP( lsr->wk[ompith][iy][ix] *=
                  lsr->kw[iy][ix] * pow(lsr->kk[iy][ix],isc); );
#else
            KOOP( lsr->wk[ompith][iy][ix] = sf_crmul(lsr->wk[ompith][iy][ix],
                                        lsr->kw[iy][ix] * pow(lsr->kk[iy][ix],isc)); );
#endif

            // FFT 
            wexfft(false,(kiss_fft_cpx**) lsr->wk[ompith],lsr->f2d[ompith]);

#ifdef SF_HAS_COMPLEX_H     
            LOOP( lsr->wk[ompith][iy][ix] *= iwdz * conjf( bw[iy][ix] ); );

            LOOP( ps[iy][ix] +=
                  lsr->wk[ompith][iy][ix] * lsr->csc[isc]; );
#else
            LOOP( lsr->wk[ompith][iy][ix] = sf_cmul(lsr->wk[ompith][iy][ix],
                                       sf_cmul(iwdz,conjf( bw[iy][ix] ))); );

            LOOP( ps[iy][ix] = sf_cadd(ps[iy][ix],
                                       sf_crmul(lsr->wk[ompith][iy][ix],lsr->csc[isc])); );
#endif
        }
    }
}

