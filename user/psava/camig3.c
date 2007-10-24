/* 3-D common-azimuth migration/modeling using extended split-step */

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

#include "camig3.h"
#include "taper3.h"
#include "slow3.h"
#include "cam3.h"

#include "slice.h"
/*^*/

#include "weutil.h"
/*^*/

#define LOOP(a) for(ihx=0;ihx<cub->ahx.n;ihx++){ \
                for(imy=0;imy<cub->amy.n;imy++){ \
                for(imx=0;imx<cub->amx.n;imx++){ {a} }}}
#define SOOP(a) for(ily=0;ily<cub->aly.n;ily++){ \
                for(ilx=0;ilx<cub->alx.n;ilx++){ {a} }}

/*------------------------------------------------------------*/
cub3d camig3_cube(bool    verb_,
		  sf_axis amx_,
		  sf_axis amy_,
		  sf_axis amz_,
		  sf_axis ahx_,
		  sf_axis alx_,
		  sf_axis aly_,
		  sf_axis aw_,
		  sf_axis ae_,
		  float   eps_,
		  int     ompnth_,
		  int     ompchunk_
    )
/*< initialize SR migration space >*/
{
    cub3d cub;
    cub = (cub3d) sf_alloc(1,sizeof(*cub));

    cub->verb=verb_;
    cub->amx = sf_nod(amx_);
    cub->amy = sf_nod(amy_);
    cub->amz = sf_nod(amz_);

    cub->ahx = sf_nod(ahx_);

    cub->alx = sf_nod(alx_);
    cub->aly = sf_nod(aly_);
    cub->ae  = sf_nod(ae_);

    cub->aw  = sf_nod(aw_);
    cub->aw.d *= 2.*SF_PI; /* from hertz to radians */
    cub->aw.o *= 2.*SF_PI;

    cub->eps      = eps_;

    cub->ompnth   = ompnth_;
    cub->ompchunk = ompchunk_;

    return cub;
}

/*------------------------------------------------------------*/
camoperator3d camig3_init(cub3d cub)
/*< initialize >*/
{
    camoperator3d weop;
    weop = (camoperator3d) sf_alloc(1,sizeof(*weop));

    weop->ww = sf_complexalloc4(cub->amx.n,cub->amy.n,cub->ahx.n,cub->ompnth);
    weop->qq = sf_floatalloc3  (cub->amx.n,cub->amy.n,cub->ahx.n);

    return weop;
}

/*------------------------------------------------------------*/
void camig3_close(camoperator3d weop)
/*< free allocated storage >*/
{
    free( ***weop->ww); free( **weop->ww);  free( *weop->ww); free( weop->ww);
    ;                   free( **weop->qq);  free( *weop->qq); free( weop->qq);
}

/*------------------------------------------------------------*/
void camig3(camoperator3d weop,
	    cub3d cub,
	    cam3d cam,
	    tap3d tap,
	    slo3d slo,
	    bool  inv   /*      forward/adjoint flag */, 
	    fslice data /* data  [nw][nhx][nmy][nmx] */,
	    fslice imag /* image [nz][nhx][nmy][nmx] */)
/*< Apply migration/modeling >*/
{
    int imz,iw,imy,imx,ihx;
    sf_complex w;
    int ompith=0;

    if (!inv) { /* prepare image for migration */
	LOOP( weop->qq[ihx][imy][imx] = 0.0; );
	for (imz=0; imz<cub->amz.n; imz++) {
	    fslice_put(imag,imz,weop->qq[0][0]);
	}
    }
    
    for (iw=0; iw<cub->aw.n; iw++) {
	if(cub->verb) sf_warning ("(ith=%d) ... <iw=%3d of %3d>",
				  ompith,iw+1,cub->aw.n);

	if (inv) { /* MODELING */
	    w = sf_cmplx(cub->eps*cub->aw.d,
			 +(cub->aw.o+iw*cub->aw.d)); /* causal */

	    LOOP( weop->ww[ompith][ihx][imy][imx] = sf_cmplx(0,0); );  
	    
	    /* upward continuation */
	    fslice_get(slo->slice,cub->amz.n-1,slo->so[ompith][0]);
	    slow3_twoway(cub,slo,slo->so,ompith); /* 2way time */

	    for (imz=cub->amz.n-1; imz>0; imz--) {

#ifdef _OPENMP	    
#pragma omp critical
#endif 
		{
		fslice_get(imag,imz,weop->qq[0][0]); /* I.C. */
#ifdef SF_HAS_COMPLEX_H
		LOOP( weop->ww[ompith][ihx][imy][imx] +=
		      weop->qq        [ihx][imy][imx]; );
#else
		LOOP( weop->ww[ompith][ihx][imy][imx].r += 
		      weop->qq        [ihx][imy][imx]; );
#endif		
		}
		fslice_get(slo->slice,imz-1,slo->ss[ompith][0]);
		slow3_twoway(cub,slo,slo->ss,ompith); /* 2way time */

		cam3_ssf(w,weop->ww[ompith],cub,cam,tap,slo,imz,ompith);	
		slow3_advance(cub,slo,ompith);
	    }

	    /* imaging at z=0 */
#ifdef _OPENMP	    
#pragma omp critical
#endif
	    {
		fslice_get(imag,0,weop->qq[0][0]);
#ifdef SF_HAS_COMPLEX_H
		LOOP( weop->ww[ompith][ihx][imy][imx]  += 
		      weop->qq        [ihx][imy][imx]; );	
#else
		LOOP( weop->ww[ompith][ihx][imy][imx].r  += 
		      weop->qq        [ihx][imy][imx]; );
#endif	  
		taper3d(weop->ww[ompith],tap);
		fslice_put(data,iw,weop->ww[ompith][0][0]);    /* output data @ iz=0 */
	    }
	    
	} else { /* MIGRATION */
	    w = sf_cmplx(cub->eps*cub->aw.d,
			 -(cub->aw.o+iw*cub->aw.d)); /* anti-causal */

	    /* imaging at z=0 */
#ifdef _OPENMP
#pragma omp critical
#endif
	    fslice_get(data,iw,weop->ww[ompith][0][0]);
	    taper3d(weop->ww[ompith],tap);
	    
#ifdef _OPENMP	 
#pragma omp critical
#endif
	    {
		fslice_get(imag,0,weop->qq[0][0]);
		LOOP(;      weop->qq        [ihx][imy][imx] += 
		     crealf(weop->ww[ompith][ihx][imy][imx] ); );
		fslice_put(imag,0,weop->qq[0][0]);
	    }
	    
	    /* downward continuation */
	    fslice_get(slo->slice,0,slo->so[ompith][0]);	
	    slow3_twoway(cub,slo,slo->so,ompith); /* 2way time */
	    for (imz=0; imz<cub->amz.n-1; imz++) {
		fslice_get(slo->slice,imz+1,slo->ss[ompith][0]);
		slow3_twoway(cub,slo,slo->ss,ompith); /* 2way time */
		
		cam3_ssf(w,weop->ww[ompith],cub,cam,tap,slo,imz,ompith);
		slow3_advance(cub,slo,ompith);

#ifdef _OPENMP	 	
#pragma omp critical
#endif
		{
		    fslice_get(imag,imz+1,weop->qq[0][0]); /* I.C. */
		    LOOP(;      weop->qq        [ihx][imy][imx] += 
			 crealf(weop->ww[ompith][ihx][imy][imx] ); );
		    fslice_put(imag,imz+1,weop->qq[0][0]);
		}
	    } /* z */
	} /* else */
    } /* w */
}

