/* 3-D SSR modeling using extended SSF */

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

#include "srmod3.h"
#include "taper3.h"
#include "slow3.h"
#include "ssr3.h"

#include "weutil.h"
/*^*/

#define LOOP(a) for(imy=0;imy<cub->amy.n;imy++){ \
                for(imx=0;imx<cub->amx.n;imx++){ \
		    {a} \
		}} /* loop in x-domain */

/*------------------------------------------------------------*/
cub3d srmod3_cube(bool    verb_,
		  sf_axis amx_,
		  sf_axis amy_,
		  sf_axis amz_,
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
ssroperator3d srmod3_init(cub3d cub)
/*< initialize SR modeling >*/
{
    ssroperator3d weop;
    weop = (ssroperator3d) sf_alloc(1,sizeof(*weop));
    
    weop->ww_s = sf_complexalloc3(cub->amx.n,cub->amy.n,cub->ompnth);
    weop->ww_r = sf_complexalloc3(cub->amx.n,cub->amy.n,cub->ompnth);

    weop->rr   = sf_floatalloc3(cub->amx.n,cub->amy.n,cub->ompnth);
    weop->wtmp = sf_fslice_init   (cub->amx.n*cub->amy.n,cub->amz.n*cub->ompnth,sizeof(sf_complex));

    return weop;
}

/*------------------------------------------------------------*/
void srmod3_close(ssroperator3d weop)
/*< free allocated storage >*/
{
    free(**weop->ww_s); free( *weop->ww_s); free( weop->ww_s);
    free(**weop->ww_r); free( *weop->ww_r); free( weop->ww_r);
    
    free( **weop->rr); free( *weop->rr); free( weop->rr);
    sf_fslice_close(weop->wtmp);
}

/*------------------------------------------------------------*/
void srmod3(ssroperator3d weop,
	    cub3d cub,
	    ssr3d ssr,
	    tap3d tap,
	    slo3d s_s,
	    slo3d s_r,
	    sf_fslice swfl /* source   data [nw][ny][nx] */,
	    sf_fslice rwfl /* receiver data [nw][ny][nx] */,
	    sf_fslice refl
    )
/*< apply SR modeling >*/
{
    int imz,iw,imx,imy,ie;
    sf_complex w;
    int ompith=0;
    int kth;
    
    for (ie=0; ie<cub->ae.n; ie++) {
	
#ifdef _OPENMP
#pragma omp parallel for schedule(static)		\
    private(ompith,iw,w,imx,imy,imz,kth)		\
    shared(swfl,rwfl,ie,weop,cub,ssr,tap,s_s,s_r)
#endif
	for (iw=0; iw<cub->aw.n; iw++) {
#ifdef _OPENMP	    
	    ompith=omp_get_thread_num();
#pragma omp critical
#endif
	    if(cub->verb) sf_warning ("(ith=%d) ... <iw=%3d of %3d> ... <ie=%3d of %3d>",
				      ompith,iw+1,cub->aw.n,ie+1,cub->ae.n);
	    
	    kth = ompith*cub->amz.n; /* tmp file thread I/O shift */
	    w = sf_cmplx(cub->eps*cub->aw.d,cub->aw.o+iw*cub->aw.d);
	    
	    /*------------------------------------------------------------*/
	    /* source wavefield */
#ifdef _OPENMP	    
#pragma omp critical
#endif    
	    sf_fslice_get(swfl,ie*cub->aw.n+iw,weop->ww_s[ompith][0]); 
	    taper2d(weop->ww_s[ompith],tap);
#ifdef _OPENMP	    
#pragma omp critical
#endif    
	    sf_fslice_put(weop->wtmp,kth+0,weop->ww_s[ompith][0]);
	    
#ifdef _OPENMP	    
#pragma omp critical
#endif	    
	    sf_fslice_get(s_s->slice,0,s_s->so[ompith][0]);
	    for (imz=0; imz<cub->amz.n-1; imz++) {
#ifdef _OPENMP	    
#pragma omp critical
#endif	    
		sf_fslice_get(s_s->slice,imz+1,s_s->ss[ompith][0]);
		
		ssr3_ssf(w,weop->ww_s[ompith],cub,ssr,tap,s_s,imz,ompith);
		
		slow3_advance(cub,s_s,ompith);
		
#ifdef _OPENMP	    
#pragma omp critical
#endif    
		sf_fslice_put(weop->wtmp,kth+imz+1,weop->ww_s[ompith][0]);
	    } /* z (down-going) */
	    
	    /*------------------------------------------------------------*/
	    /* receiver wavefield */
	    LOOP( weop->ww_r[ompith][imy][imx] = sf_cmplx(0.,0.); );

#ifdef _OPENMP	    
#pragma omp critical
#endif	    
	    sf_fslice_get(s_r->slice,cub->amz.n-1,s_r->so[ompith][0]);
	    for (imz=cub->amz.n-1; imz>0; imz--) {
#ifdef _OPENMP	    
#pragma omp critical
#endif	    
		sf_fslice_get(s_r->slice,imz-1,s_r->ss[ompith][0]);
		
#ifdef _OPENMP	    
#pragma omp critical
#endif   
		sf_fslice_get(weop->wtmp,kth+imz,weop->ww_s[ompith][0]); 
		
		sf_fslice_get(refl,imz,weop->rr[ompith][0]);
		
#ifdef SF_HAS_COMPLEX_H
		LOOP( weop->ww_s[ompith][imy][imx] *= weop->rr[ompith]  [imy][imx];
		      weop->ww_r[ompith][imy][imx] += weop->ww_s[ompith][imy][imx]; 
                    );
#else
		LOOP( weop->ww_s[ompith][imy][imx] = 
		      sf_crmul(weop->ww_s[ompith][imy][imx],
			       weop->rr[ompith][imy][imx]);
		      weop->ww_r[ompith][imy][imx] = 
		      sf_cadd(weop->ww_r[ompith][imy][imx],
			      weop->ww_s[ompith][imy][imx]); );
#endif
		
		ssr3_ssf(w,weop->ww_r[ompith],cub,ssr,tap,s_r,imz,ompith);
		
		slow3_advance(cub,s_r,ompith);
	    } /* z (up-going) */
	    
	    taper2d(weop->ww_r[ompith],tap);
	    
#ifdef _OPENMP	    
#pragma omp critical
#endif  
	    sf_fslice_put(rwfl,ie*cub->aw.n+iw,weop->ww_r[ompith][0]);
	} /* w */
	
    } /* e */
}

