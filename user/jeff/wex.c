/* 3-D SSR migration/modeling using extended split-step */

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

#ifdef _OPENMP
#include <omp.h>
#include "omputil.h"
#endif
/*^*/

#include "wex.h"
#include "wextap.h"
#include "wexslo.h"
#include "wexssr.h"
#include "wexutl.h"
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
wexcub3d wex_cube(bool    verb_,
		  sf_axis amx_,
		  sf_axis amy_,
		  sf_axis az_,
		  sf_axis alx_,
		  sf_axis aly_,
		  sf_axis aw_,
		  sf_axis ae_,
		  float   eps_,
		  int     ompnth_
    )
/*< initialize SR migration space >*/
{
    wexcub3d cub;
    cub = (wexcub3d) sf_alloc(1,sizeof(*cub));

    cub->verb=verb_;
    cub->amx = sf_nod(amx_);
    cub->amy = sf_nod(amy_);
    cub->az  = sf_nod(az_);

    cub->alx = sf_nod(alx_);
    cub->aly = sf_nod(aly_);
    cub->ae  = sf_nod(ae_);

    cub->aw    = sf_nod(aw_);
    cub->aw.d *= 2.*SF_PI; /* from hertz to radians */
    cub->aw.o *= 2.*SF_PI;

    cub->eps      = eps_;

    cub->ompnth   = ompnth_;

    return cub;
}

/*------------------------------------------------------------*/
wexop3d wex_init(wexcub3d cub)
/*< initialize >*/
{
    wexop3d weop;
    weop = (wexop3d) sf_alloc(1,sizeof(*weop));

    weop->wws = sf_complexalloc3(cub->amx.n,cub->amy.n,cub->ompnth);
    weop->wwr = sf_complexalloc3(cub->amx.n,cub->amy.n,cub->ompnth);
    weop->ws  = sf_complexalloc4(cub->amx.n,cub->amy.n,cub->az.n,cub->ompnth);
    weop->wr  = sf_complexalloc4(cub->amx.n,cub->amy.n,cub->az.n,cub->ompnth);

    weop->rr  = sf_floatalloc3(cub->amx.n,cub->amy.n,cub->az.n);

    return weop;
}


/*------------------------------------------------------------*/
void wex_store(wexop3d weop,
	       wexcub3d cub,
               int iz,
               int ompith
    )
/*< store wavefield >*/
{
    int imx,imy;
    LOOP( weop->wr [ompith][iz][imy][imx] = 
	  weop->wwr[ompith]    [imy][imx]; );
}

/*------------------------------------------------------------*/
void wex_close(wexop3d weop)
/*< free allocated storage >*/
{
                     free( **weop->wws);free(*weop->wws);free( weop->wws);
                     free( **weop->wwr);free(*weop->wwr);free( weop->wwr);
                     free( **weop->rr );free(*weop->rr );free( weop->rr );

    free(***weop->ws); free(**weop->ws); free(*weop->ws); free( weop->ws);
    free(***weop->wr); free(**weop->wr); free(*weop->wr); free( weop->wr);
}

/*------------------------------------------------------------*/
void wexdatum(wexop3d weop,
              wexcub3d cub,
              wexssr3d ssr,
              wextap3d tap,
              wexslo3d slo,
              int  wexsign,
              sf_file data /* adj=0, reflectively r[nz][nmy][nmx]; adj=1 data [nw][nmy][nmx] */,
              sf_file wflr /* receiver wavefield [nw][nz][nmy][nmx] */,
              int     inv /* down/upward flag */)
/*< Datuming wavefield by continuation >*/
{
    int iz, iw, ie;
    sf_complex w;
    int ompith=0;

    /*------------------------------------------------------------*/
    for(ie=0; ie<cub->ae.n; ie++){
#ifdef _OPENMP
#pragma omp parallel for schedule(static)  \
    private(ompith,iw,w,iz)                  \
    shared(data,wflr,weop,cub,ssr,tap,slo,inv,ie)
#endif
        for (iw=0; iw<cub->aw.n; iw++) {
#ifdef _OPENMP
            ompith = omp_get_thread_num();
#endif

#ifdef _OPENMP
#pragma omp critical
#endif
            {   /* read wavefield */
                if(cub->verb) sf_warning ("DTM... (ith=%02d) ... <iw=%3d of %3d> ... <ie=%3d of %3d>",
                                          ompith,iw+1,cub->aw.n,ie+1,cub->ae.n);
                sf_seek(data,cub->amx.n*cub->amy.n*(ie*cub->aw.n+iw)*sizeof(sf_complex),SEEK_SET);
                sf_complexread(weop->wwr[ompith][0],cub->amx.n*cub->amy.n,data);
                wextap2D(weop->wwr[ompith],tap);
            }

            if(inv) { /* UPWARD DATUMING */
                w = sf_cmplx( cub->eps*cub->aw.d, wexsign * (cub->aw.o+iw*cub->aw.d) );

                for (iz=cub->az.n-1; iz>0; iz--){
                    /* extrapolate wavefield */
                    wexssr(w,weop->wwr[ompith],cub,ssr,tap,slo,iz,ompith,false);
                }
            } else { /* DOWNWARD DATUMING */
                w = sf_cmplx( cub->eps*cub->aw.d, wexsign * (cub->aw.o+iw*cub->aw.d) );

                for (iz=0; iz<cub->az.n-1; iz++){
                    /* extrapolate wavefield */
                    wexssr(w,weop->wwr[ompith],cub,ssr,tap,slo,iz,ompith,true);
                }
            }

            wextap2D(weop->wwr[ompith],tap);

#ifdef _OPENMP      
#pragma omp critical
#endif
            {
                sf_seek(wflr,sizeof(sf_complex)*cub->amx.n*cub->amy.n*(cub->aw.n*ie+iw),SEEK_SET);
                sf_complexwrite(weop->wwr[ompith][0],cub->amx.n*cub->amy.n,wflr);
            }
        }
    }
}

/*------------------------------------------------------------*/
void wex(wexop3d weop,
	 wexcub3d cub,
	 wexssr3d ssr,
	 wextap3d tap,
	 wexslo3d slo,
	 int  wexsign,
	 sf_file   data /* adj=0, reflectively r[nz][nmy][nmx]; adj=1 data [nw][nmy][nmx] */,
	 sf_file   wflr /* receiver wavefield [nw][nz][nmy][nmx] */,
         bool    adj /* adjoint flag             */)
/*< Save wavefield from continuation >*/
{
    int imx, imy, iz, iw;
    sf_complex ****wt, ***wwt;
    sf_complex w;
    int ompith=0;

    wt = sf_complexalloc4(cub->amx.n,cub->amy.n,cub->az.n,cub->ompnth);
    wwt= sf_complexalloc3(cub->amx.n,cub->amy.n,cub->aw.n);

    if(adj)
        sf_complexread(wwt[0][0],cub->amx.n*cub->amy.n*cub->aw.n,data);
 
    /*------------------------------------------------------------*/    
#ifdef _OPENMP
#pragma omp parallel for schedule(static)  \
    private(ompith,iw,w,imx,imy,iz)			\
    shared(data,wflr,weop,cub,ssr,tap,slo,adj)
#endif
    //for (iw=0; iw<1; iw++) {
    for (iw=0; iw<cub->aw.n; iw++) {
#ifdef _OPENMP
	ompith = omp_get_thread_num();
#endif
        if(adj){ /* extrapolation for migration */
 	    
#ifdef _OPENMP
#pragma omp critical
#endif
            {
	        if(cub->verb) sf_warning ("dWR... (ith=%02d) ... <iw=%3d of %3d>",
				      ompith,iw+1,cub->aw.n);	
                LOOP( weop->wwr[ompith][imy][imx] = wwt[iw][imy][imx]; );
            }

            w = sf_cmplx( cub->eps*cub->aw.d, wexsign * (cub->aw.o+iw*cub->aw.d) );

            /* store wavefield */
            wextap2D(weop->wwr[ompith],tap);
            wex_store(weop,cub,0,ompith);

            /*------------------------------------------------------------*/
            for (iz=0; iz<cub->az.n-1; iz++) {

       //        sf_warning("bwf=%f+%f,iz=%d",crealf(weop->wwr[ompith][0][100]),cimagf(weop->wwr[ompith][0][100]),iz);

                /* extrapolate wavefield */
	        wexssr(w,weop->wwr[ompith],cub,ssr,tap,slo,iz,ompith,adj);

                /* store wavefield */
	        wextap2D(weop->wwr[ompith],tap);
	        wex_store(weop,cub,iz+1,ompith);

	    } /* z */
	    /*------------------------------------------------------------*/
	    /* output wavefield */
#ifdef _OPENMP
#pragma omp critical
#endif
            {
                sf_seek(wflr,sizeof(sf_complex)*cub->amx.n*cub->amy.n*cub->az.n*iw,SEEK_SET);
                sf_complexwrite(weop->wr[ompith][0][0],cub->amx.n*cub->amy.n*cub->az.n,wflr);
            }
        } /* end if adj */

        else{ /* extrapolation for modeling */

            /*------------------------------------------------------------*/
	    w = sf_cmplx( cub->eps*cub->aw.d, wexsign * (cub->aw.o+iw*cub->aw.d) );
            LOOP( weop->wwr[ompith][imy][imx] = sf_cmplx(0.,0.); );

#ifdef _OPENMP      
#pragma omp critical
#endif
	    {
	        if(cub->verb) sf_warning ("uWR... (ith=%02d) ... <iw=%3d of %3d>",
					  ompith,iw+1,cub->aw.n);
                sf_seek(data,cub->amx.n*cub->amy.n*cub->az.n*iw*sizeof(sf_complex),SEEK_SET);
                sf_complexread(wt[ompith][0][0],cub->amx.n*cub->amy.n*cub->az.n,data);
	    }

            for (iz=cub->az.n-1; iz>0; iz--) {

	        /* inject wavefield */ 
#ifdef SF_HAS_COMPLEX_H
	        LOOP( weop->wwr[ompith][imy][imx] += wt[ompith][iz][imy][imx]; );
#else
		LOOP( weop->wwr[ompith][imy][imx] = sf_cadd(weop->wwr[ompith][imy][imx],wt[ompith][iz][imy][imx]); );
#endif
//                sf_warning("rwf=%1.10f+%1.10f,pwf=%1.10f+%1.10f,iz=%d",crealf(weop->wwr[ompith][0][100]),cimagf(weop->wwr[ompith][0][100]),crealf(wt[ompith][iz][0][100]),cimagf(wt[ompith][iz][0][100]),iz);
                /* store wavefield */
                wextap2D(weop->wwr[ompith],tap);
                wex_store(weop,cub,iz,ompith);

		/* extrapolate wavefield */
                wexssr(w,weop->wwr[ompith],cub,ssr,tap,slo,iz,ompith,adj);

            } /* z (up-going) */
            /*------------------------------------------------------------*/
#ifdef SF_HAS_COMPLEX_H
            LOOP( weop->wwr[ompith][imy][imx] += wt[ompith][iz][imy][imx]; );
#else
            LOOP( weop->wwr[ompith][imy][imx] = sf_cadd(weop->wwr[ompith][imy][imx],wt[ompith][iz][imy][imx]); );
#endif
            wextap2D(weop->wwr[ompith],tap);
            wex_store(weop,cub,0,ompith);

            /* output wavefield */
#ifdef _OPENMP      
#pragma omp critical
#endif
            {
                sf_seek(wflr,sizeof(sf_complex)*cub->amx.n*cub->amy.n*cub->az.n*iw,SEEK_SET);
                sf_complexwrite(weop->wr[ompith][0][0],cub->amx.n*cub->amy.n*cub->az.n,wflr);
            }
        }
    } /* w */
    /*------------------------------------------------------------*/

    free( ***wt); free(  **wt);  free(   *wt);  free(    wt);
                  free(  **wwt); free(   *wwt); free(    wwt);
}


