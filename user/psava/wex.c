/* 3-D SSR extrapolation using extended split-step */

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

#ifdef _OPENMP
#include <omp.h>
#endif

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
    cub->az = sf_nod(az_);

    cub->alx = sf_nod(alx_);
    cub->aly = sf_nod(aly_);
    cub->ae  = sf_nod(ae_);

    cub->aw  = sf_nod(aw_);
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

    weop->ww = sf_complexalloc3(cub->amx.n,cub->amy.n,cub->ompnth);
    weop->w  = sf_complexalloc4(cub->amx.n,cub->amy.n,cub->az.n,cub->ompnth);

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
    LOOP( weop->w [ompith][iz][imy][imx] = 
	  weop->ww[ompith]    [imy][imx]; );
}

/*------------------------------------------------------------*/
void wex_close(wexop3d weop)
/*< free allocated storage >*/
{
                     free( **weop->ww);free(*weop->ww);free( weop->ww);
    free(***weop->w); free(**weop->w); free(*weop->w); free( weop->w);
}

/*------------------------------------------------------------*/
void wex(wexop3d weop,
	 wexcub3d cub,
	 wexssr3d ssr,
	 wextap3d tap,
	 wexslo3d slo,
	 int  wexsign,
	 sf_fslice data /*      data [nw][nmy][nmx] */,
	 sf_fslice wfld /* wavefield [nw][nmy][nmx] */)
/*< Save wavefield from downward continuation >*/
{
    int iz,iw; /* imx,imy; */
    sf_complex w;
    int ompith=0;

    /*------------------------------------------------------------*/    
#ifdef _OPENMP
#pragma omp parallel for \
    private(ompith,iw,w,iz) \
    shared(data,wfld,weop,cub,ssr,tap,slo)
#endif
    for (iw=0; iw<cub->aw.n; iw++) {
#ifdef _OPENMP
	ompith = omp_get_thread_num();
#pragma omp critical
#endif
	{
	    if(cub->verb) sf_warning ("(ith=%2d) ... <iw=%3d of %3d>",
				      ompith,iw+1,cub->aw.n);
	    
	    /* read wavefield */
	    sf_fslice_get(data,iw,weop->ww[ompith][0]);
	}
	    w = sf_cmplx( cub->eps*cub->aw.d, wexsign * (cub->aw.o+iw*cub->aw.d) );

	    /* store wavefield */
	    wextap(weop->ww[ompith],tap);
	    wex_store(weop,cub,0,ompith);

	    /*------------------------------------------------------------*/
	    for (iz=0; iz<cub->az.n-1; iz++) {

		/* extrapolate wavefield */
		wexssr(w,weop->ww[ompith],cub,ssr,tap,slo,iz,ompith);

		/* store wavefield */
		wextap(weop->ww[ompith],tap);
		wex_store(weop,cub,iz+1,ompith);

	    } /* z */
	    /*------------------------------------------------------------*/

	    /* output wavefield */
#ifdef _OPENMP
#pragma omp critical
#endif
	    sf_fslice_put(wfld,iw,weop->w[ompith][0][0]);

    } /* w */
}

/*
 
            for(imy=0; imy<cub->amy.n; imy++){
                for(imx=(cub->amx.n/2-1); imx<(cub->amx.n/2+2); imx++){
                  sf_warning("wwr=%f+%f",crealf(weop->ww[ompith][imy][imx]),cimagf(weop->ww[ompith][imy][imx])); 
                }
            }
                    
                if(iz<2){
                    for(imy=0; imy<cub->amy.n; imy++){
                        for(imx=(cub->amx.n/2-1); imx<(cub->amx.n/2+2); imx++){
                          sf_warning("iz=%d,wwr=%f+%f",iz,crealf(weop->ww[ompith][imy][imx]),cimagf(weop->ww[ompith][imy][imx])); 
                        } 
                    } 
                }  
*/

