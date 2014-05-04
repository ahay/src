/* 3-D SR MVA using extended split-step */
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
#ifdef _OPENMP
#include <omp.h>
#include "omputil.h"
#endif

#include "wex.h"
#include "wextap.h"
#include "wexeic.h"
#include "wexmva.h"
#include "wexssr.h"
#include "wexlsr.h"
#include "wexslo.h"
#include "wexutl.h"
/*^*/

#define  XOOP(a) for( iz=0;  iz<cub->az.n;   iz++) { \
                 for( imy=0; imy<cub->amy.n; imy++){ \
                 for( imx=0; imx<cub->amx.n; imx++){ \
                     {a} \
                 }}}

#define  LOOP(a) for( imy=0; imy<cub->amy.n; imy++){ \
                 for( imx=0; imx<cub->amx.n; imx++){ \
                     {a} \
                 }}

/*------------------------------------------------------------*/
wexmvaop3d wexmva_init(wexcub3d cub,
                       wexcip3d cip 
                            )
/*< initialize >*/
{

    wexmvaop3d wexmvaop;
    wexmvaop = (wexmvaop3d) sf_alloc(1,sizeof(*wexmvaop));

    /* allocate storage */
    wexmvaop->bws = sf_complexalloc4(cub->amx.n,cub->amy.n,cub->az.n,cub->ompnth);
    wexmvaop->bwr = sf_complexalloc4(cub->amx.n,cub->amy.n,cub->az.n,cub->ompnth);

    wexmvaop->dws = sf_complexalloc3(cub->amx.n,cub->amy.n,cub->ompnth);
    wexmvaop->dwr = sf_complexalloc3(cub->amx.n,cub->amy.n,cub->ompnth);
    wexmvaop->stmp= sf_complexalloc3(cub->amx.n,cub->amy.n,cub->ompnth);
    wexmvaop->rtmp= sf_complexalloc3(cub->amx.n,cub->amy.n,cub->ompnth);
    wexmvaop->pws = sf_complexalloc4(cub->amx.n,cub->amy.n,cub->az.n,cub->ompnth);
    wexmvaop->pwr = sf_complexalloc4(cub->amx.n,cub->amy.n,cub->az.n,cub->ompnth);

    wexmvaop->pss = sf_complexalloc3(cub->amx.n,cub->amy.n,cub->az.n);
    wexmvaop->psr = sf_complexalloc3(cub->amx.n,cub->amy.n,cub->az.n);

    wexmvaop->pssum = sf_complexalloc3(cub->amx.n,cub->amy.n,cub->az.n);

    return wexmvaop;
}

/*------------------------------------------------------------*/
void wexmva_close(wexmvaop3d wexmvaop)
/*< free allocated storage >*/
{
    free(**wexmvaop->bws); free( *wexmvaop->bws); free( wexmvaop->bws);
    free(**wexmvaop->bwr); free( *wexmvaop->bwr); free( wexmvaop->bwr);

    free(**wexmvaop->dws); free( *wexmvaop->dws); free( wexmvaop->dws);
    free(**wexmvaop->dwr); free( *wexmvaop->dwr); free( wexmvaop->dwr);

    free(**wexmvaop->stmp); free( *wexmvaop->stmp); free( wexmvaop->stmp);
    free(**wexmvaop->rtmp); free( *wexmvaop->rtmp); free( wexmvaop->rtmp);

    free(**wexmvaop->pss); free( *wexmvaop->pss); free( wexmvaop->pss);
    free(**wexmvaop->psr); free( *wexmvaop->psr); free( wexmvaop->psr);

    free(**wexmvaop->pssum); free( *wexmvaop->pssum); free( wexmvaop->pssum);

    free(***wexmvaop->pws); free( **wexmvaop->pws); free(  *wexmvaop->pws); free(   wexmvaop->pws);
    free(***wexmvaop->pwr); free( **wexmvaop->pwr); free(  *wexmvaop->pwr); free(   wexmvaop->pwr);

}

/*------------------------------------------------------------*/
void wexmva(wexmvaop3d      wexmvaop,
            int  adj,         /* forward/adjoint flag         */
            wexcub3d cub,     /* wavefield hypercube          */
            wexssr3d ssr,     /* SSR operator                 */
            wexlsr3d lsr,     /* SSR operator                 */
            wextap3d tap,     /* tapering operator            */
            wexslo3d slo,     /* background slowness          */
            sf_file wfls,     /* source wavefield             */
            sf_file wflr,     /* receiver wavefield           */
            sf_file pswf,   /* perturbed source wavefield   */
            sf_file prwf,   /* perturbed reciever wavefield */
            sf_file pslo      /* slowness perturbation        */
          )
/*< Forward/Adjoint WEMVA operator for non-zero offset>*/
{

    int imx, imy, iz, iw;
    int ompith = 0;
    sf_complex wr;

    if(!adj){
        sf_complexread(wexmvaop->pss[0][0],cub->amx.n*cub->amy.n*cub->az.n,pslo);
        XOOP(wexmvaop->psr[iz][imy][imx] = wexmvaop->pss[iz][imy][imx]; );
    }

#ifdef _OPENMP
#pragma omp parallel for schedule(static) \
    private(ompith,iw,wr,iz,imy,imx)       \
    shared(adj,wfls,wflr,pswf,prwf,wexmvaop,cub,ssr,lsr,tap,slo,pslo)
#endif
    /* loop over frequencies w */
    for (iw=0; iw<cub->aw.n; iw++) {
#ifdef _OPENMP
        ompith=omp_get_thread_num();
#endif

#ifdef _OPENMP
#pragma omp critical
#endif
        {
            if(cub->verb) sf_warning ("MVA (ith=%2d) ... <iw=%3d of %3d>",
                                    ompith,iw+1,cub->aw.n);
            /* load background wavefields */
            sf_seek(wfls,cub->amx.n*cub->amy.n*cub->az.n*iw*sizeof(sf_complex),SEEK_SET);
            sf_seek(wflr,cub->amx.n*cub->amy.n*cub->az.n*iw*sizeof(sf_complex),SEEK_SET);
            sf_complexread(wexmvaop->bws[ompith][0][0],cub->amx.n*cub->amy.n*cub->az.n,wfls);
            sf_complexread(wexmvaop->bwr[ompith][0][0],cub->amx.n*cub->amy.n*cub->az.n,wflr);
        }

		/* adjoint: perturbed WFs -> slowness */
        if (adj) { 
            /* ws = sf_cmplx(cub->eps*cub->aw.d,-(cub->aw.o+iw*cub->aw.d)); anti-causal */
            wr = sf_cmplx(cub->eps*cub->aw.d,+(cub->aw.o+iw*cub->aw.d)); /*      causal */

            for (iz=0; iz<=cub->az.n-1; iz++) {
                LOOP( wexmvaop->bws[ompith][iz][imy][imx] = conjf(wexmvaop->bws[ompith][iz][imy][imx]); );
            }

            LOOP( wexmvaop->dws[ompith][imy][imx]=sf_cmplx(0.,0.);
                  wexmvaop->dwr[ompith][imy][imx]=sf_cmplx(0.,0.); 
                  wexmvaop->stmp[ompith][imy][imx]=sf_cmplx(0.,0.); 
                  wexmvaop->rtmp[ompith][imy][imx]=sf_cmplx(0.,0.); );
 
#ifdef _OPENMP
#pragma omp critical
#endif
                {
                  sf_seek(pswf,cub->amx.n*cub->amy.n*cub->az.n*iw*sizeof(sf_complex),SEEK_SET);
                  sf_seek(prwf,cub->amx.n*cub->amy.n*cub->az.n*iw*sizeof(sf_complex),SEEK_SET);
                  sf_complexread(wexmvaop->pws[ompith][0][0],cub->amx.n*cub->amy.n*cub->az.n,pswf);
                  sf_complexread(wexmvaop->pwr[ompith][0][0],cub->amx.n*cub->amy.n*cub->az.n,prwf);
                }

			// Loop over Depth for both forward and backward
           for (iz=cub->az.n-1; iz>=0; iz--) {

                /* CONJ(Ts(m)) += CONJ(dWs(m)), Tr(m) += dWr(m) */
                LOOP( 
                      wexmvaop->dws[ompith][imy][imx] += wexmvaop->pws[ompith][iz][imy][imx];
                      wexmvaop->dwr[ompith][imy][imx] += wexmvaop->pwr[ompith][iz][imy][imx]; 
                    );

                /* Wavefield Scattering (W.S.) dW -> dS */                wexlsr_w2s(wr,wexmvaop->bws[ompith][iz],cub,lsr,slo,wexmvaop->dws[ompith],wexmvaop->stmp[ompith],iz,ompith);
                wexlsr_w2s(wr,wexmvaop->bwr[ompith][iz],cub,lsr,slo,wexmvaop->dwr[ompith],wexmvaop->rtmp[ompith],iz,ompith);

                /* store perturbed slowness, ds(m,z) += S^+[Ws(m,z),CONJ(Ts(m))]; ds(m) += S^+[Wr(m,z),Tr(m)] */
                LOOP(
                      wexmvaop->pssum[iz][imy][imx] += wexmvaop->stmp[ompith][imy][imx];
                      wexmvaop->pssum[iz][imy][imx] += wexmvaop->rtmp[ompith][imy][imx];
                    );

                /* wavefield extrapolation, CONJ(Ts(m)) = E^- [CONJ(Ts(m))], Tr(m) = E^+ [Tr(m)] */
                if(iz>0){
                    wexssr(wr,wexmvaop->dws[ompith],cub,ssr,tap,slo,iz,ompith,false);
                    wexssr(wr,wexmvaop->dwr[ompith],cub,ssr,tap,slo,iz,ompith,false);  
                }

            } /* loop over z */
        }
        /* forward: slowness -> perturbed WFs */
        else {   
            /* ws = sf_cmplx(cub->eps*cub->aw.d,+(cub->aw.o+iw*cub->aw.d));      causal */
            wr = sf_cmplx(cub->eps*cub->aw.d,-(cub->aw.o+iw*cub->aw.d)); /* anti-causal */

            for (iz=0; iz<=cub->az.n-1; iz++) {
                LOOP( wexmvaop->bws[ompith][iz][imy][imx] = conjf(wexmvaop->bws[ompith][iz][imy][imx]); );
            }
       
      		XOOP( wexmvaop->pws[ompith][iz][imy][imx]=sf_cmplx(0.,0.);
                  wexmvaop->pwr[ompith][iz][imy][imx]=sf_cmplx(0.,0.); );

            LOOP( wexmvaop->dws[ompith][imy][imx]=sf_cmplx(0.,0.);
                  wexmvaop->dwr[ompith][imy][imx]=sf_cmplx(0.,0.); 
                );

            for (iz=0; iz<=cub->az.n-1; iz++) {
            
                /* Wavefield Scattering (W.S.) dS -> dW */
                wexlsr_s2w(wr,wexmvaop->bws[ompith][iz],cub,lsr,slo,wexmvaop->stmp[ompith],wexmvaop->pss[iz],iz,ompith);
                wexlsr_s2w(wr,wexmvaop->bwr[ompith][iz],cub,lsr,slo,wexmvaop->rtmp[ompith],wexmvaop->psr[iz],iz,ompith);

                /* CONJ(Ts(m)) += S^- [CONJ(Ws(m)),ds(m,z)], Tr(m) += S^- [Wr(m),ds(m,z)] */
                LOOP( wexmvaop->dws[ompith][imy][imx] += wexmvaop->stmp[ompith][imy][imx]; );
                LOOP( wexmvaop->dwr[ompith][imy][imx] += wexmvaop->rtmp[ompith][imy][imx]; );

                /* store perturbed wavefield, dWs(m,z) = Ts(m), dWr(m,z) = Tr(m) */
                LOOP( 
                      wexmvaop->pws[ompith][iz][imy][imx] = wexmvaop->dws[ompith][imy][imx];
                      wexmvaop->pwr[ompith][iz][imy][imx] = wexmvaop->dwr[ompith][imy][imx];
                    );   
                    
                /* wavefield extrapolation, Ts(m) = E+ [CONJ(Ts(m))], Tr(m) = E- [Tr(m)] */
                if(iz<cub->az.n-1){
                    wexssr(wr,wexmvaop->dws[ompith],cub,ssr,tap,slo,iz,ompith,true);
                    wexssr(wr,wexmvaop->dwr[ompith],cub,ssr,tap,slo,iz,ompith,true);
                }

            } /* loop over z */
#ifdef _OPENMP
#pragma omp critical
#endif
            { /* put perturbed wavefields */
                sf_seek(pswf,sizeof(sf_complex)*cub->amx.n*cub->amy.n*cub->az.n*iw,SEEK_SET);
                sf_seek(prwf,sizeof(sf_complex)*cub->amx.n*cub->amy.n*cub->az.n*iw,SEEK_SET);
                sf_complexwrite(wexmvaop->pws[ompith][0][0],(size_t)cub->amx.n*cub->amy.n*cub->az.n,pswf);
                sf_complexwrite(wexmvaop->pwr[ompith][0][0],(size_t)cub->amx.n*cub->amy.n*cub->az.n,prwf);
            }
        } /* loop over if adj else */
    } /* loop over w */
         
    if(adj)
        sf_complexwrite(wexmvaop->pssum[0][0],cub->amx.n*cub->amy.n*cub->az.n,pslo);
}

/*------------------------------------------------------------*/
void wexzomva(wexmvaop3d      wexmvaop,
            int  adj,         /* forward/adjoint flag         */
            wexcub3d cub,     /* wavefield hypercube          */
            wexssr3d ssr,     /* SSR operator                 */
            wexlsr3d lsr,     /* SSR operator                 */
            wextap3d tap,     /* tapering operator            */
            wexslo3d slo,     /* background slowness          */
            sf_file wflr,     /* receiver wavefield           */
            sf_file prwf,   /* perturbed reciever wavefield */
            sf_file pslo      /* slowness perturbation        */
          )
/*< Forward/Adjoint WEMVA operator for zero offset >*/
{

    int imx, imy, iz, iw;
    int ompith = 0;
    sf_complex wr;

    if(!adj){
        sf_complexread(wexmvaop->psr[0][0],cub->amx.n*cub->amy.n*cub->az.n,pslo);
    }

#ifdef _OPENMP
#pragma omp parallel for schedule(static) \
    private(ompith,iw,wr,iz,imy,imx)       \
    shared(adj,wflr,prwf,wexmvaop,cub,ssr,lsr,tap,slo,pslo)
#endif
    /* loop over frequencies w */
    for (iw=0; iw<cub->aw.n; iw++) {
#ifdef _OPENMP
        ompith=omp_get_thread_num();
#endif

#ifdef _OPENMP
#pragma omp critical
#endif
        {
            if(cub->verb) sf_warning ("MVA (ith=%2d) ... <iw=%3d of %3d>",
                                    ompith,iw+1,cub->aw.n);
            /* load background wavefields */
            sf_seek(wflr,cub->amx.n*cub->amy.n*cub->az.n*iw*sizeof(sf_complex),SEEK_SET);
            sf_complexread(wexmvaop->bwr[ompith][0][0],cub->amx.n*cub->amy.n*cub->az.n,wflr);
        }

        if (adj) { /* adjoint: perturbed WFs -> slowness */
            wr = sf_cmplx(cub->eps*cub->aw.d,+(cub->aw.o+iw*cub->aw.d)); /*      causal */

            LOOP( wexmvaop->dwr[ompith][imy][imx]=sf_cmplx(0.,0.); );
 
#ifdef _OPENMP
#pragma omp critical
#endif
                {
                  sf_seek(prwf,cub->amx.n*cub->amy.n*cub->az.n*iw*sizeof(sf_complex),SEEK_SET);
                  sf_complexread(wexmvaop->pwr[ompith][0][0],cub->amx.n*cub->amy.n*cub->az.n,prwf);
                }

           for (iz=cub->az.n-1; iz>=0; iz--) {
                /* CONJ(Ts(m)) += CONJ(dWs(m)), Tr(m) += dWr(m) */
                LOOP( wexmvaop->dwr[ompith][imy][imx] += wexmvaop->pwr[ompith][iz][imy][imx];);

                /* Wavefield Scattering (W.S.) dW -> dS */
          		wexlsr_w2s(wr,wexmvaop->bwr[ompith][iz],cub,lsr,slo,wexmvaop->dwr[ompith],wexmvaop->rtmp[ompith],iz,ompith);

                /* store perturbed slowness, ds(m,z) += S^+[Ws(m,z),CONJ(Ts(m))]; ds(m) += S^+[Wr(m,z),Tr(m)] */
                LOOP( wexmvaop->pssum[iz][imy][imx] += wexmvaop->rtmp[ompith][imy][imx]; );

                /* wavefield extrapolation, CONJ(Ts(m)) = E^- [CONJ(Ts(m))], Tr(m) = E^+ [Tr(m)] */
                if(iz>0){
                    wexssr(wr,wexmvaop->dwr[ompith],cub,ssr,tap,slo,iz,ompith,false);  
                }

            } /* loop over z */
        }
        else {   /* forward: slowness -> perturbed WFs */
            wr = sf_cmplx(cub->eps*cub->aw.d,-(cub->aw.o+iw*cub->aw.d)); /* anti-causal */
//
            XOOP( wexmvaop->pwr[ompith][iz][imy][imx]=sf_cmplx(0.,0.); );
            LOOP( wexmvaop->dwr[ompith][imy][imx]=sf_cmplx(0.,0.);     );

            for (iz=0; iz<=cub->az.n-1; iz++) {
                /* Wavefield Scattering (W.S.) dS -> dW */
                wexlsr_s2w(wr,wexmvaop->bwr[ompith][iz],cub,lsr,slo,wexmvaop->rtmp[ompith],wexmvaop->psr[iz],iz,ompith);

                /* CONJ(Ts(m)) += S^- [CONJ(Ws(m)),ds(m,z)], Tr(m) += S^- [Wr(m),ds(m,z)] */
                LOOP( wexmvaop->dwr[ompith][imy][imx] += wexmvaop->rtmp[ompith][imy][imx]; );

                /* store perturbed wavefield, dWs(m,z) = Ts(m), dWr(m,z) = Tr(m) */
                LOOP( wexmvaop->pwr[ompith][iz][imy][imx] = wexmvaop->dwr[ompith][imy][imx];);   
                /* wavefield extrapolation, Ts(m) = E+ [CONJ(Ts(m))], Tr(m) = E- [Tr(m)] */
                if(iz<cub->az.n-1){
                    wexssr(wr,wexmvaop->dwr[ompith],cub,ssr,tap,slo,iz,ompith,true);
                }

            } /* loop over z */
#ifdef _OPENMP
#pragma omp critical
#endif
            { /* put perturbed wavefields */
                sf_seek(prwf,sizeof(sf_complex)*cub->amx.n*cub->amy.n*cub->az.n*iw,SEEK_SET);
                sf_complexwrite(wexmvaop->pwr[ompith][0][0],(size_t)cub->amx.n*cub->amy.n*cub->az.n,prwf);
            }
        } /* loop over if adj else */
    } /* loop over w */
         
    if(adj)
        sf_complexwrite(wexmvaop->pssum[0][0],cub->amx.n*cub->amy.n*cub->az.n,pslo);
}
