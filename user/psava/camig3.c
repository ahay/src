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

#ifdef _OPENMP
#include <omp.h>
#endif

#include "camig3.h"
#include "taper3.h"
#include "slow3.h"
#include "cam3.h"

#include "weutil.h"
/*^*/

#define LOOP(a) for(ihx=0;ihx<cub->ahx.n;ihx++){ \
                for(imy=0;imy<cub->amy.n;imy++){ \
                for(imx=0;imx<cub->amx.n;imx++){ \
		    {a} \
		}}} /* loop in x-domain */

/*------------------------------------------------------------*/
cub3d camig3_cube(bool   verb_,
		  sf_axis amx_,
		  sf_axis amy_,
		  sf_axis amz_,
		  sf_axis ahx_,
		  sf_axis alx_,
		  sf_axis aly_,
		  sf_axis aw_,
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
	    sf_fslice data /* data  [nw][nhx][nmy][nmx] */,
	    sf_fslice imag /* image [nz][nhx][nmy][nmx] */)
/*< Apply migration/modeling >*/
{
    int imz,iw,imy,imx,ihx;
    sf_complex w;
    int ompith=0;

    if (!inv) { /* prepare image for migration */
	LOOP( weop->qq[ihx][imy][imx] = 0.0; );
	for (imz=0; imz<cub->amz.n; imz++) {
	    sf_fslice_put(imag,imz,weop->qq[0][0]);
	}
    }
    
#ifdef _OPENMP
#pragma omp parallel for schedule(static)	\
    private(ompith,iw,w,imx,imy,imz,ihx)	\
    shared(data,weop,cub,cam,tap,slo)
#endif
    for (iw=0; iw<cub->aw.n; iw++) {
#ifdef _OPENMP
	ompith=omp_get_thread_num();
#endif

	if (inv) { /* MODELING */
	    w = sf_cmplx(cub->eps*cub->aw.d,+(cub->aw.o+iw*cub->aw.d)); /* causal */

	    LOOP( weop->ww[ompith][ihx][imy][imx] = sf_cmplx(0,0); );  
	    
	    /* upward continuation */
#ifdef _OPENMP
#pragma omp critical
#endif
	    sf_fslice_get(slo->slice,cub->amz.n-1,slo->so[ompith][0]);
	    for (imz=cub->amz.n-1; imz>0; imz--) {
#ifdef _OPENMP
#pragma omp critical
#endif
		if(cub->verb) sf_warning ("(ith=%d) ... <iw=%3d of %3d> <iz=%3d of %3d>",
					  ompith,iw+1,cub->aw.n,imz+1,cub->amz.n);

#ifdef _OPENMP	    
#pragma omp critical
#endif 
		{
		sf_fslice_get(imag,imz,weop->qq[0][0]); /* I.C. */
#ifdef SF_HAS_COMPLEX_H
		LOOP( weop->ww[ompith][ihx][imy][imx] +=
		      weop->qq        [ihx][imy][imx]; );
#else
		LOOP( weop->ww[ompith][ihx][imy][imx].r += 
		      weop->qq        [ihx][imy][imx]; );
#endif		
		}
#ifdef _OPENMP
#pragma omp critical
#endif
		sf_fslice_get(slo->slice,imz-1,slo->ss[ompith][0]);
		cam3_ssf(w,weop->ww[ompith],cub,cam,tap,slo,imz,ompith);	
		slow3_advance(cub,slo,ompith);
	    }

	    /* imaging at z=0 */
#ifdef _OPENMP	    
#pragma omp critical
#endif
	    {
		sf_fslice_get(imag,0,weop->qq[0][0]);
#ifdef SF_HAS_COMPLEX_H
		LOOP( weop->ww[ompith][ihx][imy][imx]  += 
		      weop->qq        [ihx][imy][imx]; );	
#else
		LOOP( weop->ww[ompith][ihx][imy][imx].r  += 
		      weop->qq        [ihx][imy][imx]; );
#endif	  
		taper3d(weop->ww[ompith],tap);
		sf_fslice_put(data,iw,weop->ww[ompith][0][0]);    /* output data @ iz=0 */
	    }
	    
	} else { /* MIGRATION */
	    w = sf_cmplx(cub->eps*cub->aw.d,-(cub->aw.o+iw*cub->aw.d)); /* anti-causal */

	    /* imaging at z=0 */
#ifdef _OPENMP
#pragma omp critical
#endif
	    sf_fslice_get(data,iw,weop->ww[ompith][0][0]);
	    taper3d(weop->ww[ompith],tap);
	    
#ifdef _OPENMP	 
#pragma omp critical
#endif
	    {
		sf_fslice_get(imag,0,weop->qq[0][0]);
		LOOP(;      weop->qq        [ihx][imy][imx] += 
		     crealf(weop->ww[ompith][ihx][imy][imx] ); );
		sf_fslice_put(imag,0,weop->qq[0][0]);
	    }
	    
	    /* downward continuation */
#ifdef _OPENMP
#pragma omp critical
#endif
	    sf_fslice_get(slo->slice,0,slo->so[ompith][0]);	
	    for (imz=0; imz<cub->amz.n-1; imz++) {
#ifdef _OPENMP
#pragma omp critical
#endif
		if(cub->verb) sf_warning ("(ith=%d) ... <iw=%3d of %3d> <iz=%3d of %3d>",
					  ompith,iw+1,cub->aw.n,imz+1,cub->amz.n);
#ifdef _OPENMP
#pragma omp critical
#endif
		sf_fslice_get(slo->slice,imz+1,slo->ss[ompith][0]);
		cam3_ssf(w,weop->ww[ompith],cub,cam,tap,slo,imz,ompith);
		slow3_advance(cub,slo,ompith);

#ifdef _OPENMP	 	
#pragma omp critical
#endif
		{
		    sf_fslice_get(imag,imz+1,weop->qq[0][0]); /* I.C. */
		    LOOP(;      weop->qq        [ihx][imy][imx] += 
			 crealf(weop->ww[ompith][ihx][imy][imx] ); );
		    sf_fslice_put(imag,imz+1,weop->qq[0][0]);
		}

	    } /* z */

	} /* else */
    } /* w */
}

/*------------------------------------------------------------*/
void cadtm3(camoperator3d weop,
	    cub3d cub,
	    cam3d cam,
	    tap3d tap,
	    slo3d slo,
	    bool  inv   /* forward/adjoint flag */, 
	    sf_fslice data /* data [nw][nmy][nmx] */,
	    sf_fslice wfld /* wfld [nw][nmy][nmx] */)
/*< Apply upward/downward datuming >*/
{
    int imz,iw;
    sf_complex w;
    int ompith=0;
    
#ifdef _OPENMP
#pragma omp parallel for schedule(static)	\
    private(ompith,iw,w,imz)			\
    shared(data,wfld,weop,cub,cam,tap,slo)
#endif
    for (iw=0; iw<cub->aw.n; iw++) {
#ifdef _OPENMP
	ompith=omp_get_thread_num();
#endif
	
	if(inv) { /* UPWARD DATUMING */
	    w = sf_cmplx(cub->eps*cub->aw.d,+(cub->aw.o+iw*cub->aw.d)); /* causal */

#ifdef _OPENMP	
#pragma omp critical
#endif
	    sf_fslice_get(wfld,iw,weop->ww[ompith][0][0]);
	    taper3d(weop->ww[ompith],tap);
	
#ifdef _OPENMP
#pragma omp critical
#endif
	    sf_fslice_get(slo->slice,cub->amz.n-1,slo->so[ompith][0]);
	    for (imz=cub->amz.n-1; imz>0; imz--) {
#ifdef _OPENMP
#pragma omp critical
#endif
		if(cub->verb) sf_warning ("(ith=%d) ... <iw=%3d of %3d> <iz=%3d of %3d>",
					  ompith,iw+1,cub->aw.n,imz+1,cub->amz.n);

#ifdef _OPENMP
#pragma omp critical
#endif
		sf_fslice_get(slo->slice,imz-1,slo->ss[ompith][0]);		
		cam3_ssf(w,weop->ww[ompith],cub,cam,tap,slo,imz,ompith);
		slow3_advance(cub,slo,ompith);
	    }
	    taper3d(weop->ww[ompith],tap);
#ifdef _OPENMP	
#pragma omp critical
#endif
	    sf_fslice_put(data,iw,weop->ww[ompith][0][0]);

	} else { /* DOWNWARD DATUMING */
	    w = sf_cmplx(cub->eps*cub->aw.d,-(cub->aw.o+iw*cub->aw.d)); /* anti-causal */
	    
#ifdef _OPENMP	
#pragma omp critical
#endif
	    sf_fslice_get(data,iw,weop->ww[ompith][0][0]);
	    taper3d(weop->ww[ompith],tap);
	    
#ifdef _OPENMP
#pragma omp critical
#endif
	    sf_fslice_get(slo->slice,0,slo->so[ompith][0]);
	    for (imz=0; imz<cub->amz.n-1; imz++) {
#ifdef _OPENMP
#pragma omp critical
#endif
		if(cub->verb) sf_warning ("(ith=%d) ... <iw=%3d of %3d> <iz=%3d of %3d>",
					  ompith,iw+1,cub->aw.n,imz+1,cub->amz.n);

#ifdef _OPENMP
#pragma omp critical
#endif
		sf_fslice_get(slo->slice,imz+1,slo->ss[ompith][0]);
		cam3_ssf(w,weop->ww[ompith],cub,cam,tap,slo,imz,ompith);
		slow3_advance(cub,slo,ompith);
	    }
	    taper3d(weop->ww[ompith],tap);
#ifdef _OPENMP	
#pragma omp critical
#endif
	    sf_fslice_put(wfld,iw,weop->ww[ompith][0][0]);
	} /* else */
    } /* w */
    
}

/*------------------------------------------------------------*/
void cawfl3(camoperator3d weop,
	    cub3d cub,
	    cam3d cam,
	    tap3d tap,
	    slo3d slo,
	    bool  inv   /* forward/adjoint flag */, 
	    sf_fslice data /*      data [nw][nmy][nmx] */,
	    sf_fslice wfld /* wavefield [nw][nmy][nmx] */)
/*< Save wavefield from downward continuation >*/
{
    int imz,iw;
    sf_complex w;
    int ompith=0;
    
#ifdef _OPENMP
#pragma omp parallel for schedule(static)	\
    private(ompith,iw,w,imz)			\
    shared(data,wfld,weop,cub,cam,tap,slo)
#endif
    for (iw=0; iw<cub->aw.n; iw++) {
#ifdef _OPENMP
	ompith=omp_get_thread_num();
#endif

	if(inv) { /*   UPWARD EXTRAPOLATION */
	    w = sf_cmplx(cub->eps*cub->aw.d,+(cub->aw.o+iw*cub->aw.d)); /* causal */
	    
#ifdef _OPENMP	
#pragma omp critical
#endif
	    sf_fslice_get(data,iw,weop->ww[ompith][0][0]);
	    taper3d(weop->ww[ompith],tap);
	    sf_fslice_put(wfld,iw*cub->amz.n+cub->amz.n-1,weop->ww[ompith][0][0]);

#ifdef _OPENMP
#pragma omp critical
#endif
	    sf_fslice_get(slo->slice,cub->amz.n-1,slo->so[ompith][0]);
	    for (imz=cub->amz.n-1; imz>0; imz--) {
#ifdef _OPENMP
#pragma omp critical
#endif
		if(cub->verb) sf_warning ("(ith=%d) ... <iw=%3d of %3d> <iz=%3d of %3d>",
					  ompith,iw+1,cub->aw.n,imz+1,cub->amz.n);

#ifdef _OPENMP
#pragma omp critical
#endif
		sf_fslice_get(slo->slice,imz-1,slo->ss[ompith][0]);
		cam3_ssf(w,weop->ww[ompith],cub,cam,tap,slo,imz,ompith);
		slow3_advance(cub,slo,ompith);

#ifdef _OPENMP	
#pragma omp critical
#endif
		sf_fslice_put(wfld,iw*cub->amz.n+imz-1,weop->ww[ompith][0][0]);
	    }

	} else {  /* DOWNWARD EXTRAPOLATION */
	    w = sf_cmplx(cub->eps*cub->aw.d,-(cub->aw.o+iw*cub->aw.d)); /* anti-causal */

#ifdef _OPENMP	
#pragma omp critical
#endif
	    sf_fslice_get(data,iw,weop->ww[ompith][0][0]);
	    taper3d(weop->ww[ompith],tap);
	    sf_fslice_put(wfld,iw*cub->amz.n,weop->ww[ompith][0][0]);

#ifdef _OPENMP
#pragma omp critical
#endif
	    sf_fslice_get(slo->slice,0,slo->so[ompith][0]);
	    for (imz=0; imz<cub->amz.n-1; imz++) {	
#ifdef _OPENMP
#pragma omp critical
#endif
		if(cub->verb) sf_warning ("(ith=%d) ... <iw=%3d of %3d> <iz=%3d of %3d>",
					  ompith,iw+1,cub->aw.n,imz+1,cub->amz.n);
    
#ifdef _OPENMP
#pragma omp critical
#endif
		sf_fslice_get(slo->slice,imz+1,slo->ss[ompith][0]);
		cam3_ssf(w,weop->ww[ompith],cub,cam,tap,slo,imz,ompith);
		slow3_advance(cub,slo,ompith);

#ifdef _OPENMP	
#pragma omp critical
#endif
		sf_fslice_put(wfld,iw*cub->amz.n+imz+1,weop->ww[ompith][0][0]);
	    } /* z */
	} /* else */
    } /* w */
}
