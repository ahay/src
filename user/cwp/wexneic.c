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
#include "wexneic.h"
/*^*/
#ifdef SF_HAS_COMPLEX_H
#define cWGH(a,b,c)  (+1.*(conjf(a)*b)*c)
#define cAWGH(a,b,c) (+1.*(a*b)*c)
#define cCOR(a,b) (conjf(a)*b)
#define cACOR(a,b) (a*b)
#define cMUL(a,b) (a*b) 
#else
#define cWGH(a,b,c) (+1.*sf_cmul((sf_cmul(conjf(a),b)),c))
#define cAWGH(a,b,c) (+1.*sf_cmul((sf_cmul(a,b)),c))
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
wexncip3d wexncip_init(wexcub3d cub,
                        int  nhx_,
                        int  nhy_,
                        int  nhz_,
                        int  nht_,
                        int  nhx2_,
                        int  nhy2_,
                        int  nhz2_,
                        int  nht2_,
                        int   nc_,
		       float  dhx_,
		       float  dhy_,
		       float  dhz_,
		       float  dht_,
		       float  oht_,
		       sf_file Fc,
		       sf_file Fdx,
		       sf_file Fdy,
		       int eic)
/*< initialize I.C. >*/
{

    int iw, ic, ihx, ihy, ihz, iht, icx, icy, icz;
    float ht, w, cosAx, cosAy, sinAx, sinAy;

    float px, py, pxz, pyz, nx, ny, nz;
    int pxs, pxzs, pys, pyzs, nxs, nys, nzs;
    /*------------------------------------------------------------*/
    wexncip3d cip;
    cip = (wexncip3d) sf_alloc(1,sizeof(*cip));

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
	cip->dhx = dhx_;
	cip->dhy = dhy_;
	cip->dhz = dhz_;

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

        for(ic=0; ic<cip->nc; ic++){
            for(iht=0; iht<cip->nht2; iht++){
                for(ihz=0; ihz<cip->nhz2; ihz++){
                    for(ihy=0; ihy<cip->nhy2; ihy++){
                        for(ihx=0; ihx<cip->nhx2; ihx++){
                            cip->ei[ic][iht][ihz][ihy][ihx] = sf_cmplx(0.0,0.0);
                        }
                    }
                }
            }
        }

        /*------------------------------------------------------------*/
        /* CIP coordinates */
        cip->cc= (pt3d*) sf_alloc(cip->nc,sizeof(*cip->cc));
        pt3dread1(Fc,cip->cc,cip->nc,3); /* read coordinates */

        /* dip field */
        cip->dipx = sf_floatalloc3(cub->amx.n,cub->amy.n,cub->az.n);
        sf_floatread(cip->dipx[0][0],cub->amx.n*cub->amy.n*cub->az.n,Fdx);

	cip->dipy = sf_floatalloc3(cub->amx.n,cub->amy.n,cub->az.n);
        if(cub->amy.n>1){
	  sf_floatread(cip->dipy[0][0],cub->amx.n*cub->amy.n*cub->az.n,Fdy);
	}

        cip->mcx=sf_intalloc4(cip->nhx2,cip->nhy2,cip->nhz2,cip->nc);
        cip->pcx=sf_intalloc4(cip->nhx2,cip->nhy2,cip->nhz2,cip->nc);
        cip->mcy=sf_intalloc4(cip->nhx2,cip->nhy2,cip->nhz2,cip->nc);
        cip->pcy=sf_intalloc4(cip->nhx2,cip->nhy2,cip->nhz2,cip->nc);
        cip->mcz=sf_intalloc4(cip->nhx2,cip->nhy2,cip->nhz2,cip->nc);
        cip->pcz=sf_intalloc4(cip->nhx2,cip->nhy2,cip->nhz2,cip->nc);

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
	      icy = 0.5+(cip->cc[ic].y-cub->amy.o)/cub->amy.d;
	      icx = 0.5+(cip->cc[ic].x-cub->amx.o)/cub->amx.d;
	      icz = 0.5+(cip->cc[ic].z-cub->az.o)/cub->az.d;
	  
              cosAx = 1.0/sqrt(1.0+cip->dipx[icz][icy][icx]*cip->dipx[icz][icy][icx]);
              if(cip->dipx[icz][icy][icx]<0)
                sinAx = -1.0*sqrt(1.0-cosAx*cosAx);
              else
                sinAx =  1.0*sqrt(1.0-cosAx*cosAx);

              if(cub->amy.n>1){
                cosAy = 1.0/sqrt(1.0+cip->dipy[icz][icy][icx]*cip->dipy[icz][icy][icx]);
                if(cip->dipy[icz][icy][icx]<0)
                  sinAy = -1.0*sqrt(1.0-cosAy*cosAy);
                else
                  sinAy =  1.0*sqrt(1.0-cosAy*cosAy);
              }
              else{
                cosAy = 1.0;
                sinAy = 0.0;
              }

              sf_warning("cosAx = %f, sinAx = %f, cosAy = %f, sinAy = %f",cosAx, sinAx, cosAy, sinAy);

	      for(ihy=-cip->nhy; ihy<cip->nhy+1; ihy++){
		py  = cip->dhy*ihy*cosAy;
		pys = py/cub->amy.d;
		
		pyz  = cip->dhy*ihy*sinAy;
		pyzs = pyz/cub->az.d;
		
	        for(ihx=-cip->nhx; ihx<cip->nhx+1; ihx++){
		  px  = cip->dhx*ihx*cosAx;
		  pxs = px/cub->amx.d;

		  pxz  = cip->dhx*ihx*sinAx;
		  pxzs = pxz/cub->az.d;

		  for(ihz=-cip->nhz; ihz<cip->nhz+1; ihz++){
		    nx  = cip->dhz*ihz*(1.0*sinAx*cosAy);
		    nxs = nx/cub->amx.d;
		  
		    ny  = cip->dhz*ihz*(1.0*cosAx*sinAy);
                    nys = ny/cub->amy.d;
                   
		    nz = cip->dhz*ihz*(-1.0*cosAx*cosAy);
		    nzs = nz/cub->az.d;

		    cip->mcx[ic][cip->nhz+ihz][cip->nhy+ihy][cip->nhx+ihx] = icx+nxs+pxs;
		    cip->pcx[ic][cip->nhz+ihz][cip->nhy+ihy][cip->nhx+ihx] = icx-nxs-pxs;

		    cip->mcy[ic][cip->nhz+ihz][cip->nhy+ihy][cip->nhx+ihx] = icy+nys+pys;
		    cip->pcy[ic][cip->nhz+ihz][cip->nhy+ihy][cip->nhx+ihx] = icy-nys-pys;		  

		    cip->mcz[ic][cip->nhz+ihz][cip->nhy+ihy][cip->nhx+ihx] = icz-nzs-pxzs-pyzs;
		    cip->pcz[ic][cip->nhz+ihz][cip->nhy+ihy][cip->nhx+ihx] = icz+nzs+pxzs+pyzs;

//                    sf_warning("ihx=%d,ihy=%d,ihz=%d,pxs=%d,nxs=%d,pys=%d,nys=%d",ihx,ihy,ihz,pxs,nxs,pys,nys);
                    //sf_warning("ihx=%d,ihy=%d,ihz=%d,mcx=%d,pcx=%d,mcy=%d,pcy=%d,mcz=%d,pcz=%d",ihx,ihy,ihz,cip->mcx[ic][cip->nhz+ihz][cip->nhy+ihy][cip->nhx+ihx],cip->pcx[ic][cip->nhz+ihz][cip->nhy+ihy][cip->nhx+ihx],cip->mcy[ic][cip->nhz+ihz][cip->nhy+ihy][cip->nhx+ihx],cip->pcy[ic][cip->nhz+ihz][cip->nhy+ihy][cip->nhx+ihx],cip->mcz[ic][cip->nhz+ihz][cip->nhy+ihy][cip->nhx+ihx],cip->pcz[ic][cip->nhz+ihz][cip->nhy+ihy][cip->nhx+ihx]);
		  }
		}
	      }
	    }
	} /* loop over nc */  
    }

    /*for(ihx=0; ihx<cip->nhx2; ihx++){
      sf_warning("mcx=%d, mcz=%d, pcx=%d, pcz=%d",cip->mcx[0][cip->nhz][ihx],cip->mcz[0][cip->nhz][ihx],cip->pcx[0][cip->nhz][ihx],cip->pcz[0][cip->nhz][ihx]);
      }*/

    return cip;
}


/*------------------------------------------------------------*/
void wexncip_close(wexncip3d cip, int eic)
/*< free allocated storage >*/
{

    if(eic){
        free(****cip->ei); free( ***cip->ei); free(  **cip->ei); free(   *cip->ei); free(    cip->ei);
        free(*cip->tt); free( cip->tt);

        free(***cip->mcx); free(**cip->mcx), free(*cip->mcx); free(cip->mcx);
        free(***cip->pcx); free(**cip->pcx), free(*cip->pcx); free(cip->pcx);
        free(***cip->mcy); free(**cip->mcy), free(*cip->mcy); free(cip->mcy);
        free(***cip->pcy); free(**cip->pcy), free(*cip->pcy); free(cip->pcy);   
        free(***cip->mcz); free(**cip->mcz), free(*cip->mcz); free(cip->mcz);
        free(***cip->pcz); free(**cip->pcz), free(*cip->pcz); free(cip->pcz);

        free(cip->ccin);
        free(cip->cc);
    }

    free( **cip->ci);  free(  *cip->ci);  free(   cip->ci);
    free( **cip->dipx); free(  *cip->dipx); free(   cip->dipx);
    free( **cip->dipy); free(  *cip->dipy); free(   cip->dipy);

    free( **cip->ws); free(  *cip->ws); free(   cip->ws);
    free( **cip->wr); free(  *cip->wr); free(   cip->wr);

}

/*------------------------------------------------------------*/
void wexncip_for(wexcub3d  cub,
                 wexncip3d cip,
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
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) \
    private(ic,ihx,ihy,ihz,iht,wt,mcx,mcy,mcz,pcx,pcy,pcz,uu) \
    shared(cip)
#endif
	  for(ic=0; ic<cip->nc; ic++) { // loop over CIPs 
	    if(cip->ccin[ic]) {
	      for      (ihx=0; ihx<cip->nhx2; ihx++) { 		
		for    (ihy=0; ihy<cip->nhy2; ihy++) { 
		  for  (ihz=0; ihz<cip->nhz2; ihz++) { 
		    mcx = cip->mcx[ic][ihz][ihy][ihx]; pcx = cip->pcx[ic][ihz][ihy][ihx];
		    mcy = cip->mcy[ic][ihz][ihy][ihx]; pcy = cip->pcy[ic][ihz][ihy][ihx];
		    mcz = cip->mcz[ic][ihz][ihy][ihx]; pcz = cip->pcz[ic][ihz][ihy][ihx];

   		    uu = cCOR(cip->ws[mcz][mcy][mcx],cip->wr[pcz][pcy][pcx]);
		    //sf_warning("ihx=%d, ihz=%d, mcx=%d, pcx=%d, mcz=%d, pcz=%d, ws=%f, wr=%f,r=%f",ihx, ihz, mcx,pcx,mcz,pcz,crealf(cip->ws[mcz][mcy][mcx]),crealf(cip->wr[pcz][pcy][pcx]),crealf(uu));
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
        } /* end if eic */
        else{/* zero-lag imaging condition */
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) \
    private(iz,imy,imx) \
    shared(cip,cub)
#endif
            for(     iz=0;   iz<cub->az.n;  iz++){
              for(  imy=0; imy<cub->amy.n; imy++){  
                for(imx=0; imx<cub->amx.n; imx++){
                  cip->ci[iz][imy][imx] += cCOR(cip->ws[iz][imy][imx],cip->wr[iz][imy][imx]); 
                } /* loop over amx */
              } /* loop over amy */
            } /* loop over az */
        } /* end else eic */
    } /* loop over w */

    if(eic)
        sf_complexwrite(cip->ei[0][0][0][0],cip->nhx2*cip->nhy2*cip->nhz2*cip->nht2*cip->nc,imag);
    else
        sf_complexwrite(cip->ci[0][0],cub->amx.n*cub->amy.n*cub->az.n,imag);
}

//		  sf_warning("ihx=%d, ihz=%d, mcx=%d, pcx=%d, mcz=%d, pcz=%d",ihx, ihz, cip->mcx[ic][cip->nhz+ihz][cip->nhx+ihx],cip->pcx[ic][cip->nhz+ihz][cip->nhx+ihx],cip->mcz[ic][cip->nhz+ihz][cip->nhx+ihx],cip->pcz[ic][cip->nhz+ihz][cip->nhx+ihx]);
