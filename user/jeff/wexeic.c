/* OpenMP lagged-products in the time- or freq-domain */

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

/* 
 * inputs:      wavefields organized as z-x-w
 * output: lagged products organized as hz-hx-ht @ irregular points
 */

#include <rsf.h>
#include <math.h>
/*^*/
#ifdef _OPENMP
#include <omp.h>
#include "omputil.h"
#endif
/*^*/
#include "wex.h"
#include "wexeic.h"
/*^*/
#ifdef SF_HAS_COMPLEX_H
#define cWGH(a,b,c)  (+1.*(conjf(a)*b)*c)
#define cAWGH(a,b,c) (+1.*(a*b)*c)
#define cCOR(a,b) (conjf(a)*b)
#define cACOR(a,b) (a*b)
#define cMUL(a,b) (a*b) 
#else
#define cWGH(a,b,c) (sf_cmul((sf_cmul(conjf(a),b)),c))
#define cAWGH(a,b,c) (sf_cmul((sf_cmul(a,b)),c))
#define cCOR(a,b) (sf_cmul(conjf(a),b))
#define cACOR(a,b) (sf_cmul(a,b))
#define cMUL(a,b) (sf_cmul(a,b))
#endif
/*^*/
#define LOOP(a) for(imy=0;imy<cub->amy.n;imy++){ \
                for(imx=0;imx<cub->amx.n;imx++){ \
                    {a} \
                }} /* loop in x-domain */

#define  XOOP(a) for( iz=0;  iz<cub->az.n;    iz++){ \
                 for( imy=0; imy<cub->amy.n; imy++){ \
                 for( imx=0; imx<cub->amx.n; imx++){ \
                     {a} \
                 }}}

/*------------------------------------------------------------*/
wexcip3d wexcip_init(wexcub3d cub,
                        int  nhx_,
                        int  nhy_,
                        int  nhz_,
                        int  nht_,
                        int  nhx2_,
                        int  nhy2_,
                        int  nhz2_,
                        int  nht2_,
                        int   nc_,
                      float  dht_,
                      float  oht_,
                    sf_file Fc,
                        int eic)
/*< initialize I.C. >*/
{

    int iw, ic, ihx, ihy, ihz, iht, icx, icy, icz;
    float ht, w;

    /*------------------------------------------------------------*/
    wexcip3d cip;
    cip = (wexcip3d) sf_alloc(1,sizeof(*cip));

    /*------------------------------------------------------------*/
    /* allocate wavefields storage */
    cip->ws = sf_complexalloc3(cub->amx.n,cub->amy.n,cub->az.n);
    cip->wr = sf_complexalloc3(cub->amx.n,cub->amy.n,cub->az.n);
    cip->ci = sf_complexalloc3(cub->amx.n,cub->amy.n,cub->az.n);

    if(eic){
        cip->nhx = nhx_;
        cip->nhy = nhy_;
        cip->nhz = nhz_;
        cip->nht = nht_;
        cip->nhx2= nhx2_;
        cip->nhy2= nhy2_;
        cip->nhz2= nhz2_;
        cip->nht2= nht2_;
        cip->nc = nc_;

        cip->oht = oht_;
        cip->dht = dht_;

        /*------------------------------------------------------------*/
        /* precompute phase for time delay */
        cip->tt = sf_complexalloc2(cip->nht2,cub->aw.n);

        for(iw=0; iw<cub->aw.n; iw++) {
            w = -( cub->aw.o+iw*cub->aw.d );

            for(iht=0; iht<cip->nht2; iht++) {
                ht = cip->oht + (iht+0)*cip->dht;

                cip->tt[iw][iht] = sf_cmplx(cosf(2*w*ht),sinf(2*w*ht));
            }
        }

        /*------------------------------------------------------------*/
        /* allocate image storage */
        cip->ei = sf_complexalloc5(cip->nhx2,cip->nhy2,cip->nhz2,cip->nht2,cip->nc);
        cip->di = sf_complexalloc5(cip->nhx2,cip->nhy2,cip->nhz2,cip->nht2,cip->nc);

        for(ic=0; ic<cip->nc; ic++){
            for(iht=0; iht<cip->nht2; iht++){
                for(ihz=0; ihz<cip->nhz2; ihz++){
                    for(ihy=0; ihy<cip->nhy2; ihy++){
                        for(ihx=0; ihx<cip->nhx2; ihx++){
                            cip->ei[ic][iht][ihz][ihy][ihx] = sf_cmplx(0.0,0.0);
                            cip->di[ic][iht][ihz][ihy][ihx] = sf_cmplx(0.0,0.0);
                        }
                    }
                }
            }
        }

        /*------------------------------------------------------------*/
        /* CIP coordinates */
        cip->cc= (pt3d*) sf_alloc(cip->nc,sizeof(*cip->cc));
        pt3dread1(Fc,cip->cc,cip->nc,3); /* read coordinates */

        cip->mcxall=sf_intalloc2(cip->nhx2,cip->nc);
        cip->pcxall=sf_intalloc2(cip->nhx2,cip->nc);
        cip->mcyall=sf_intalloc2(cip->nhy2,cip->nc);
        cip->pcyall=sf_intalloc2(cip->nhy2,cip->nc);
        cip->mczall=sf_intalloc2(cip->nhz2,cip->nc);
        cip->pczall=sf_intalloc2(cip->nhz2,cip->nc);

        cip->ccin=sf_intalloc(cip->nc);

        cip->cxmin = cub->amx.o +               cip->nhx *cub->amx.d;
        cip->cxmax = cub->amx.o + (cub->amx.n-1-cip->nhx)*cub->amx.d;
        cip->cymin = cub->amy.o +               cip->nhy *cub->amy.d;
        cip->cymax = cub->amy.o + (cub->amy.n-1-cip->nhy)*cub->amy.d;
        cip->czmin = cub->az.o  +               cip->nhz *cub->az.d;
        cip->czmax = cub->az.o  + (cub->az.n -1-cip->nhz)*cub->az.d;

        for(ic=0; ic<cip->nc; ic++) {
            cip->ccin[ic]=(cip->cc[ic].x>=cip->cxmin && cip->cc[ic].x<=cip->cxmax &&
                           cip->cc[ic].y>=cip->cymin && cip->cc[ic].y<=cip->cymax &&
                           cip->cc[ic].z>=cip->czmin && cip->cc[ic].z<=cip->czmax)?1:0;

            if(cip->ccin[ic]) {

                icx = 0.5+(cip->cc[ic].x-cub->amx.o)/cub->amx.d;
                for(ihx=-cip->nhx; ihx<cip->nhx+1; ihx++) {
                    cip->mcxall[ic][cip->nhx+ihx] = icx-ihx;
                    cip->pcxall[ic][cip->nhx+ihx] = icx+ihx;
                }
 
                icy = 0.5+(cip->cc[ic].y-cub->amy.o)/cub->amy.d;
                for(ihy=-cip->nhy; ihy<cip->nhy+1; ihy++) {
                    cip->mcyall[ic][cip->nhy+ihy] = icy-ihy;
                    cip->pcyall[ic][cip->nhy+ihy] = icy+ihy;
                }

                icz = 0.5+(cip->cc[ic].z-cub->az.o)/cub->az.d;
                for(ihz=-cip->nhz; ihz<cip->nhz+1; ihz++) {
                    cip->mczall[ic][cip->nhz+ihz] = icz-ihz;
                    cip->pczall[ic][cip->nhz+ihz] = icz+ihz;
                }

                for(ihx=-cip->nhx; ihx<cip->nhx+1; ihx++) {
                  for(ihy=-cip->nhy; ihy<cip->nhy+1; ihy++) {
                    for(ihz=-cip->nhz; ihz<cip->nhz+1; ihz++) {
//                      sf_warning("ihx=%d,ihy=%d,ihz=%d,mcx=%d,pcx=%d,mcy=%d,pcy=%d,mcz=%d,pcz=%d",ihx,ihy,ihz,cip->mcxall[ic][cip->nhx+ihx],cip->pcxall[ic][cip->nhx+ihx],cip->mcyall[ic][cip->nhy+ihy],cip->pcyall[ic][cip->nhy+ihy],cip->mczall[ic][cip->nhz+ihz],cip->pczall[ic][cip->nhz+ihz]);
                    }

                  }
                }
            }
        } /* loop over nc */
    }

    return cip;
}


/*------------------------------------------------------------*/
void wexcip_close(wexcip3d cip, int eic)
/*< free allocated storage >*/
{

    if(eic){
        free(****cip->ei); free( ***cip->ei); free(  **cip->ei); free(   *cip->ei); free(    cip->ei);
        free(****cip->di); free( ***cip->di); free(  **cip->di); free(   *cip->di); free(    cip->di);
        free(*cip->tt); free( cip->tt);

        free(*cip->mcxall); free( cip->mcxall);
        free(*cip->mcyall); free( cip->mcyall);
        free(*cip->mczall); free( cip->mczall);
 
        free(*cip->pczall); free( cip->pczall);
        free(*cip->pcxall); free( cip->pcxall);
        free(*cip->pcyall); free( cip->pcyall);

        free(cip->ccin);
        free(cip->cc);
    }

                      free( **cip->ci); free(  *cip->ci); free(   cip->ci);

    free( **cip->ws); free(  *cip->ws); free(   cip->ws);
    free( **cip->wr); free(  *cip->wr); free(   cip->wr);

}

/*------------------------------------------------------------*/
void wexcip_for(wexcub3d cub,
                wexcip3d cip,
                 sf_file wfls,
                 sf_file wflr,
                 sf_file imag,
                 int eic,
                 int rep,
                 int cjg
         )
/*< E.I.C >*/
{

    int ihx, ihy, ihz, iht, ic, iw, iz, imy, imx;
    int mcz, mcy, mcx;
    int pcz, pcy, pcx;

    sf_complex wt, uu;
    if(rep){
        if(eic)
            sf_complexread(cip->ei[0][0][0][0],cip->nhx2*cip->nhy2*cip->nhz2*cip->nht2*cip->nc,imag);
        else
            sf_complexread(cip->ci[0][0],cub->amx.n*cub->amy.n*cub->az.n,imag);
    }

    for (iw=0; iw<cub->aw.n; iw++) {
        if(cub->verb){
          if(eic)
              sf_warning("eIC ... <iw=%3d of %3d>",iw+1,cub->aw.n);
          else
              sf_warning("cIC ... <iw=%3d of %3d>",iw+1,cub->aw.n);
        }
 
        sf_seek(wfls,cub->amx.n*cub->amy.n*cub->az.n*iw*sizeof(sf_complex),SEEK_SET);
        sf_complexread(cip->ws[0][0],cub->amx.n*cub->amy.n*cub->az.n,wfls);
        sf_seek(wflr,cub->amx.n*cub->amy.n*cub->az.n*iw*sizeof(sf_complex),SEEK_SET);
        sf_complexread(cip->wr[0][0],cub->amx.n*cub->amy.n*cub->az.n,wflr);

        if(eic){ /* extended imaging condition */
          if(cjg){
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) \
    private(ic,ihx,ihy,ihz,iht,wt,mcx,mcy,mcz,pcx,pcy,pcz,uu) \
    shared(cip)
#endif
            for(ic=0; ic<cip->nc; ic++) { // loop over CIPs 
                if(cip->ccin[ic]) {
                    for            (ihx=0; ihx<cip->nhx2; ihx++) { mcx=cip->mcxall[ic][ihx]; pcx=cip->pcxall[ic][ihx];
                        for        (ihy=0; ihy<cip->nhy2; ihy++) { mcy=cip->mcyall[ic][ihy]; pcy=cip->pcyall[ic][ihy];
                            for    (ihz=0; ihz<cip->nhz2; ihz++) { mcz=cip->mczall[ic][ihz]; pcz=cip->pczall[ic][ihz];
                                uu = cACOR(cip->ws[mcz][mcy][mcx],cip->wr[pcz][pcy][pcx]);
                                for(iht=0; iht<cip->nht2; iht++) {
                                    wt = cip->tt[iw][iht];
#ifdef SF_HAS_COMPLEX_H
         cip->ei[ic][iht][ihz][ihy][ihx] += cACOR(uu,wt); //uu);
#else
	 cip->ei[ic][iht][ihz][ihy][ihx] = sf_cadd(cip->ei[ic][iht][ihz][ihy][ihx],cACOR(uu,wt));
#endif
                                } /* loop over ht */
                            } /* loop over hz */
                        } /* loop over hy */
                    } /* loop over hx */
                } /* end if cc */
            } /* loop over nc */
          } /* end if cjg */
          else{
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) \
    private(ic,ihx,ihy,ihz,iht,wt,mcx,mcy,mcz,pcx,pcy,pcz,uu) \
    shared(cip)
#endif
            for(ic=0; ic<cip->nc; ic++) { // loop over CIPs 
                if(cip->ccin[ic]) {
                    for            (ihx=0; ihx<cip->nhx2; ihx++) { mcx=cip->mcxall[ic][ihx]; pcx=cip->pcxall[ic][ihx];
                        for        (ihy=0; ihy<cip->nhy2; ihy++) { mcy=cip->mcyall[ic][ihy]; pcy=cip->pcyall[ic][ihy];
                            for    (ihz=0; ihz<cip->nhz2; ihz++) { mcz=cip->mczall[ic][ihz]; pcz=cip->pczall[ic][ihz];
                                uu = cCOR(cip->ws[mcz][mcy][mcx],cip->wr[pcz][pcy][pcx]);
                                for(iht=0; iht<cip->nht2; iht++) {
                                    wt = cip->tt[iw][iht];
#ifdef SF_HAS_COMPLEX_H
                                    cip->ei[ic][iht][ihz][ihy][ihx] += cACOR(uu,wt);
#else
				    cip->ei[ic][iht][ihz][ihy][ihx] = sf_cadd(cip->ei[ic][iht][ihz][ihy][ihx],cACOR(uu,wt));
#endif
                                } /* loop over ht */
                            } /* loop over hz */
                        } /* loop over hy */
                    } /* loop over hx */
                } /* end if cc */
            } /* loop over nc */
          } /* end else cjg */
        } /* end if eic */
        else{/* zero-lag imaging condition */
            if(cjg){
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) \
    private(iz,imy,imx) \
    shared(cip,cub)
#endif
              for(     iz=0;   iz<cub->az.n;  iz++){
                for(  imy=0; imy<cub->amy.n; imy++){  
                  for(imx=0; imx<cub->amx.n; imx++){
#ifdef SF_HAS_COMPLEX_H
                    cip->ci[iz][imy][imx] += cACOR(cip->ws[iz][imy][imx],cip->wr[iz][imy][imx]); 
#else
		    cip->ci[iz][imy][imx] = sf_cadd(cip->ci[iz][imy][imx],cACOR(cip->ws[iz][imy][imx],cip->wr[iz][imy][imx])); 
#endif
                  } /* loop over amx */
                } /* loop over amy */
              } /* loop over az */
            } /* end if cjg */ 
            else{
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) \
    private(iz,imy,imx) \
    shared(cip,cub)
#endif
              for(     iz=0;   iz<cub->az.n;  iz++){
                for(  imy=0; imy<cub->amy.n; imy++){  
                  for(imx=0; imx<cub->amx.n; imx++){
#ifdef SF_HAS_COMPLEX_H
                    cip->ci[iz][imy][imx] += cCOR(cip->ws[iz][imy][imx],cip->wr[iz][imy][imx]); 
#else
		    cip->ci[iz][imy][imx] = sf_cadd(cip->ci[iz][imy][imx],cCOR(cip->ws[iz][imy][imx],cip->wr[iz][imy][imx])); 
#endif
                  } /* loop over ax  */
                } /* loop over amy  */
              } /* loop over az   */
            } /* end else cjg */
        } /* end else eic */
    } /* loop over w */

    if(eic)
        sf_complexwrite(cip->ei[0][0][0][0],cip->nhx2*cip->nhy2*cip->nhz2*cip->nht2*cip->nc,imag);
    else
        sf_complexwrite(cip->ci[0][0],cub->amx.n*cub->amy.n*cub->az.n,imag);
}

/*------------------------------------------------------------*/
void wexcip_adj(wexcub3d cub,
                wexcip3d cip,
                 sf_file wfls,
                 sf_file wflr,
                 sf_file imag,
                 int eic,
                 int rep,
                 int cjg
         )
/*< E.I.C >*/
{

    int ihx, ihy, ihz, iht, ic, iw, iz, imy, imx;
    int mcz, mcy, mcx;
    int pcz, pcy, pcx;
    sf_complex wt, uu;

    if(rep){
        if(eic)
            sf_complexread(cip->ei[0][0][0][0],cip->nhx2*cip->nhy2*cip->nhz2*cip->nht2*cip->nc,imag);
        else
            sf_complexread(cip->ci[0][0],cub->amx.n*cub->amy.n*cub->az.n,imag);
    }

    for (iw=0; iw<cub->aw.n; iw++) {
        XOOP( cip->wr[iz][imy][imx] = sf_cmplx(0.0,0.0); );

        if(eic)
            sf_warning("eIC ... <iw=%3d of %3d>",iw+1,cub->aw.n);
        else
            sf_warning("cIC ... <iw=%3d of %3d>",iw+1,cub->aw.n);

        sf_seek(wfls,cub->amx.n*cub->amy.n*cub->az.n*iw*sizeof(sf_complex),SEEK_SET);
        sf_complexread(cip->ws[0][0],cub->amx.n*cub->amy.n*cub->az.n,wfls);

        if(eic){
          if(cjg){
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) \
    private(ic,ihx,ihy,ihz,iht,wt,mcx,mcy,mcz,pcx,pcy,pcz,uu) \
    shared(cip)
#endif
            for(ic=0; ic<cip->nc; ic++) { // loop over CIPs
                if(cip->ccin[ic]) {
                    for            (ihx=0; ihx<cip->nhx2; ihx++) { mcx=cip->mcxall[ic][ihx]; pcx=cip->pcxall[ic][ihx];
                        for        (ihy=0; ihy<cip->nhy2; ihy++) { mcy=cip->mcyall[ic][ihy]; pcy=cip->pcyall[ic][ihy];
                            for    (ihz=0; ihz<cip->nhz2; ihz++) { mcz=cip->mczall[ic][ihz]; pcz=cip->pczall[ic][ihz];
                                uu = sf_cmplx(0.0,0.0);
                                for(iht=0; iht<cip->nht2; iht++) {
                                    wt = cip->tt[iw][iht];
#ifdef SF_HAS_COMPLEX_H
                                    cip->wr[mcz][mcy][mcx] += cWGH(cip->ws[pcz][pcy][pcx],cip->ei[ic][iht][ihz][ihy][ihx],conjf(wt));
#else
				    cip->wr[mcz][mcy][mcx] = sf_cadd(cip->wr[mcz][mcy][mcx],cWGH(cip->ws[pcz][pcy][pcx],cip->ei[ic][iht][ihz][ihy][ihx],conjf(wt)));
#endif
                  //                  sf_warning("ps=%f+%f,br=%f+%f,dr=%f+%f,iz=%d",crealf(cip->wr[mcz][mcy][mcx]),cimagf(cip->wr[mcz][mcy][mcx]),crealf(cip->ws[pcz][pcy][pcx]),cimagf(cip->ws[pcz][pcy][pcx]),crealf(cip->ei[ic][iht][ihz][ihy][ihx]),cimagf(cip->ei[ic][iht][ihz][ihy][ihx]),mcz);
                                } /* loop over ht */

                            } /* loop over hz */
                        } /* loop over hy */
                    } /* loop over hx */
                } /* end if cc */
            } /* loop over nc */
          } /* end if cjg */
          else{
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) \
    private(ic,ihx,ihy,ihz,iht,wt,mcx,mcy,mcz,pcx,pcy,pcz,uu) \
    shared(cip)
#endif
            for(ic=0; ic<cip->nc; ic++) { // loop over CIPs
                if(cip->ccin[ic]) {
                    for            (ihx=0; ihx<cip->nhx2; ihx++) { mcx=cip->mcxall[ic][ihx]; pcx=cip->pcxall[ic][ihx];
                        for        (ihy=0; ihy<cip->nhy2; ihy++) { mcy=cip->mcyall[ic][ihy]; pcy=cip->pcyall[ic][ihy];
                            for    (ihz=0; ihz<cip->nhz2; ihz++) { mcz=cip->mczall[ic][ihz]; pcz=cip->pczall[ic][ihz];
                                uu = sf_cmplx(0.0,0.0);
                                for(iht=0; iht<cip->nht2; iht++) {
                                    wt = cip->tt[iw][iht];
#ifdef SF_HAS_COMPLEX_H
                                    uu += cACOR(cip->ei[ic][iht][ihz][ihy][ihx],conjf(wt));
#else
				    uu = sf_cadd(uu,cACOR(cip->ei[ic][iht][ihz][ihy][ihx],conjf(wt)));
#endif
                                } /* loop over ht */
#ifdef SF_HAS_COMPLEX_H
                                cip->wr[pcz][pcy][pcx] += cACOR(cip->ws[mcz][mcy][mcx],uu);
#else
				cip->wr[pcz][pcy][pcx] = sf_cadd(cip->wr[pcz][pcy][pcx],cACOR(cip->ws[mcz][mcy][mcx],uu));
#endif
                            } /* loop over hz */
                        } /* loop over hy */
                    } /* loop over hx */
                } /* end if cc */
            } /* loop over nc */
          } /* end else cjg */
        } /* end if eic */

        else{
            if(cjg){
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) \
    private(iz,imy,imx) \
    shared(cip,cub)
#endif
              for(     iz=0;   iz<cub->az.n;  iz++){
                for(  imy=0; imy<cub->amy.n; imy++){  
                  for(imx=0; imx<cub->amx.n; imx++){
                    cip->wr[iz][imy][imx] = cCOR(cip->ws[iz][imy][imx],cip->ci[iz][imy][imx]);
                  } /* loop over amx */
                } /* loop over amy */
              } /* loop over az */
            } /* end if cjg */ 
            else{
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) \
    private(iz,imy,imx) \
    shared(cip,cub)
#endif
              for(     iz=0;   iz<cub->az.n;  iz++){
                for(  imy=0; imy<cub->amy.n; imy++){  
                  for(imx=0; imx<cub->amx.n; imx++){
                    cip->wr[iz][imy][imx] = cACOR(cip->ws[iz][imy][imx],cip->ci[iz][imy][imx]);
                  } /* loop over amx */
                } /* loop over amy */
              } /* loop over az */
            } /* end if cjg */ 
        } /* end if eic */

        sf_seek(wflr,sizeof(sf_complex)*cub->amx.n*cub->amy.n*cub->az.n*iw,SEEK_SET);
        sf_complexwrite(cip->wr[0][0],cub->amx.n*cub->amy.n*cub->az.n,wflr);
    } /* loop over w */
}

/*------------------------------------------------------------*/
void wexcip_for_drv(wexcub3d cub,
                   wexcip3d cip,
                    sf_file wfls,
                    sf_file wflr,
                    sf_file derv
         )
/*< E.I.C >*/
{

    int ihx, ihy, ihz, iht, ic, iw;
    int   mcz, mcy, mcx;
    int   pcz, pcy, pcx;
    sf_complex wt, dr, uu;
    float w;

    for (iw=0; iw<cub->aw.n; iw++) {

        sf_warning("eIC derivative ... <iw=%3d of %3d>",iw+1,cub->aw.n);
 
        sf_seek(wfls,cub->amx.n*cub->amy.n*cub->az.n*iw*sizeof(sf_complex),SEEK_SET);
        sf_complexread(cip->ws[0][0],cub->amx.n*cub->amy.n*cub->az.n,wfls);

        sf_seek(wflr,cub->amx.n*cub->amy.n*cub->az.n*iw*sizeof(sf_complex),SEEK_SET);
        sf_complexread(cip->wr[0][0],cub->amx.n*cub->amy.n*cub->az.n,wflr);

        w = cub->aw.o+iw*cub->aw.d;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) \
    private(ic,ihx,ihy,ihz,iht,wt,mcx,mcy,mcz,pcx,pcy,pcz,uu,dr) \
    shared(cip,w)
#endif
        for(ic=0; ic<cip->nc; ic++) { // loop over CIPs 
            if(cip->ccin[ic]) {
                for            (ihx=0; ihx<cip->nhx2; ihx++) { mcx=cip->mcxall[ic][ihx]; pcx=cip->pcxall[ic][ihx];
                    for        (ihy=0; ihy<cip->nhy2; ihy++) { mcy=cip->mcyall[ic][ihy]; pcy=cip->pcyall[ic][ihy];
                        for    (ihz=0; ihz<cip->nhz2; ihz++) { mcz=cip->mczall[ic][ihz]; pcz=cip->pczall[ic][ihz];
                            uu = cCOR(cip->ws[mcz][mcy][mcx],cip->wr[pcz][pcy][pcx]);
                            for(iht=0; iht<cip->nht2; iht++) {
                                wt = cip->tt[iw][iht];
                                dr = sf_cmplx(cub->eps*cub->aw.d,-2.0*w);
#ifdef SF_HAS_COMPLEX_H
                                dr = wt*dr;
                                cip->ei[ic][iht][ihz][ihy][ihx] += cACOR(uu,dr);
#else
                                dr = sf_cmul(wt,dr);
                                cip->ei[ic][iht][ihz][ihy][ihx] = sf_cadd(cip->ei[ic][iht][ihz][ihy][ihx],cACOR(uu,dr));
#endif

                            } /* loop over ht */
                        } /* loop over hz */
                    } /* loop over hy */
                } /* loop over hx */
            } /* end if cc */
        } /* loop over nc */
    }

    sf_complexwrite(cip->ei[0][0][0][0],cip->nhx2*cip->nhy2*cip->nhz2*cip->nht2*cip->nc,derv);
}

/*------------------------------------------------------------*/
void wexcip_for_new(wexcub3d cub,
                   wexcip3d cip,
                    sf_file wfls,
                    sf_file wflr,
                    sf_file imag
         )
/*< E.I.C >*/
{

    int ihx, ihy, ihz, iht, ic, iw;
    int   mcz, mcy, mcx;
    int   pcz, pcy, pcx;
    sf_complex wt, uu;

    for (iw=0; iw<cub->aw.n; iw++) {

        sf_warning("eIC ... <iw=%3d of %3d>",iw+1,cub->aw.n);
 
        sf_seek(wfls,cub->amx.n*cub->amy.n*cub->az.n*iw*sizeof(sf_complex),SEEK_SET);
        sf_complexread(cip->ws[0][0],cub->amx.n*cub->amy.n*cub->az.n,wfls);

        sf_seek(wflr,cub->amx.n*cub->amy.n*cub->az.n*iw*sizeof(sf_complex),SEEK_SET);
        sf_complexread(cip->wr[0][0],cub->amx.n*cub->amy.n*cub->az.n,wflr);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) \
    private(ic,ihx,ihy,ihz,iht,wt,mcx,mcy,mcz,pcx,pcy,pcz,uu) \
    shared(cip)
#endif
        for(ic=0; ic<cip->nc; ic++) { // loop over CIPs 
            if(cip->ccin[ic]) {
                for            (ihx=0; ihx<cip->nhx2; ihx++) { mcx=cip->mcxall[ic][ihx]; pcx=cip->pcxall[ic][ihx];
                    for        (ihy=0; ihy<cip->nhy2; ihy++) { mcy=cip->mcyall[ic][ihy]; pcy=cip->pcyall[ic][ihy];
                        for    (ihz=0; ihz<cip->nhz2; ihz++) { mcz=cip->mczall[ic][ihz]; pcz=cip->pczall[ic][ihz];
                            uu = cCOR(cip->ws[mcz][mcy][mcx],cip->wr[pcz][pcy][pcx]);
                            for(iht=0; iht<cip->nht2; iht++) {
                                wt = cip->tt[iw][iht];
#ifdef SF_HAS_COMPLEX_H
                                wt *= sf_cmplx(0.0,-1.0);
                                cip->ei[ic][iht][ihz][ihy][ihx] += cACOR(uu,wt);
#else
                                wt = sf_cmul(wt,sf_cmplx(0.0,-1.0));
				cip->ei[ic][iht][ihz][ihy][ihx] = sf_cadd(cip->ei[ic][iht][ihz][ihy][ihx],cACOR(uu,wt));
#endif

                            } /* loop over ht */
                        } /* loop over hz */
                    } /* loop over hy */
                } /* loop over hx */
            } /* end if cc */
        } /* loop over nc */
    } /* loop over w */

    sf_complexwrite(cip->ei[0][0][0][0],cip->nhx2*cip->nhy2*cip->nhz2*cip->nht2*cip->nc,imag);
}

/*------------------------------------------------------------*/
void wexzocip_for(wexcub3d cub,
                  wexcip3d cip,
                   sf_file wflr,
                   sf_file imag
         )
/*< zero-offset I.C >*/
{

    int iw, iz, imy, imx;

    for (iw=0; iw<cub->aw.n; iw++) {
        if(cub->verb)
            sf_warning("cIC ... <iw=%3d of %3d>",iw+1,cub->aw.n);
 
        sf_seek(wflr,cub->amx.n*cub->amy.n*cub->az.n*iw*sizeof(sf_complex),SEEK_SET);
        sf_complexread(cip->wr[0][0],cub->amx.n*cub->amy.n*cub->az.n,wflr);

        /* zero-lag imaging condition */	
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) \
    private(iz,imy,imx) \
    shared(cip,cub)
#endif
	  for(     iz=0;   iz<cub->az.n;  iz++){
	    for(  imy=0; imy<cub->amy.n; imy++){  
	      for(imx=0; imx<cub->amx.n; imx++){
#ifdef SF_HAS_COMPLEX_H
		cip->ci[iz][imy][imx] += cip->wr[iz][imy][imx]; 
#else
		cip->ci[iz][imy][imx] = sf_cadd(cip->ci[iz][imy][imx],cip->wr[iz][imy][imx]); 
#endif
	      } /* loop over amx */
	    } /* loop over amy */
	  } /* loop over az */
    } /* loop over w */
    
    sf_complexwrite(cip->ci[0][0],cub->amx.n*cub->amy.n*cub->az.n,imag);
}

/*------------------------------------------------------------*/
void wexzocip_adj(wexcub3d cub,
                wexcip3d cip,
                 sf_file wflr,
                 sf_file imag
         )
/*< zero-offset I.C >*/
{

    int iw, iz, imy, imx;

    sf_complexread(cip->ci[0][0],cub->amx.n*cub->amy.n*cub->az.n,imag);

    for (iw=0; iw<cub->aw.n; iw++) {
        XOOP( cip->wr[iz][imy][imx] = sf_cmplx(0.0,0.0); );

	if(cub->verb)
            sf_warning("cIC ... <iw=%3d of %3d>",iw+1,cub->aw.n);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) \
    private(iz,imy,imx) \
    shared(cip,cub)
#endif
	for(     iz=0;   iz<cub->az.n;  iz++){
	  for(  imy=0; imy<cub->amy.n; imy++){  
	    for(imx=0; imx<cub->amx.n; imx++){
	      cip->wr[iz][imy][imx] = cip->ci[iz][imy][imx];
	    } /* loop over amx */
	  } /* loop over amy */
	} /* loop over az */
	
        sf_seek(wflr,sizeof(sf_complex)*cub->amx.n*cub->amy.n*cub->az.n*iw,SEEK_SET);
        sf_complexwrite(cip->wr[0][0],cub->amx.n*cub->amy.n*cub->az.n,wflr);
    } /* loop over w */
}

/*


*/
