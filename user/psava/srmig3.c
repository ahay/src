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
#endif

#include "srmig3.h"
#include "taper3.h"
#include "slow3.h"
#include "ssr3.h"
#include "img3.h"

#include "slice.h"
/*^*/

#include "weutil.h"
/*^*/

/*------------------------------------------------------------*/
sroperator3d srmig3_init(bool    verb_,
			 float   eps_,
			 float   dtmax,
			 sf_axis ae_     /* experiments (e.g. shots) */,
			 sf_axis aw_     /* frequency */,
			 sf_axis amx_    /* i-line (data) */,
			 sf_axis amy_    /* x-line (data) */,
			 sf_axis amz_    /* depth */,
			 sf_axis alx_    /* i-line (slowness/image) */,
			 sf_axis aly_    /* x-line (slowness/image) */,
			 int     tmx, 
			 int     tmy     /* taper size */,
			 int     pmx, 
			 int     pmy,    /* padding in the k domain */
			 int    nrmax,   /* maximum number of references */
			 fslice slow_s_,
			 fslice slow_r_,
			 int     ompnth_
    )
/*< initialize SR migration >*/
{
    float dsmax;

    sroperator3d srop;
    srop = (sroperator3d) sf_alloc(1,sizeof(*srop));
    
    srop->verb=verb_;
    srop->eps = eps_;
    srop->ompnth = ompnth_;

    srop->amx = sf_nod(amx_);
    srop->amy = sf_nod(amy_);
    srop->amz = sf_nod(amz_);
    srop->alx = sf_nod(alx_);
    srop->aly = sf_nod(aly_);
    srop->ae  = sf_nod(ae_);

    /* from hertz to radian */
    srop->aw = sf_nod(aw_);
    srop->aw.d *= 2.*SF_PI; 
    srop->aw.o *= 2.*SF_PI;

    dsmax = dtmax/srop->amz.d;

    /* SSR */
    srop->ssr = ssr3_init(srop->amz,
			  srop->amx,
			  srop->amy,
			  srop->alx,
			  srop->aly,
			  pmx ,pmy,
			  tmx ,tmy,
			  dsmax,
			  srop->ompnth);
    
    /* taper */
    srop->tap = taper_init(srop->amx.n,srop->amy.n,1,
			   SF_MIN(tmx,srop->amx.n-1),SF_MIN(tmy,srop->amy.n-1),0,
			   true,true,false);

    /* slowness */
    srop->s_s = slow3_init(slow_s_,srop->alx,srop->aly,srop->amz,nrmax,dsmax,srop->ompnth);
    srop->s_r = slow3_init(slow_r_,srop->alx,srop->aly,srop->amz,nrmax,dsmax,srop->ompnth);
    
    /* allocate wavefield storage */
    srop->ww_s = sf_complexalloc3(srop->amx.n,srop->amy.n,srop->ompnth);
    srop->ww_r = sf_complexalloc3(srop->amx.n,srop->amy.n,srop->ompnth);

    return srop;
}

/*------------------------------------------------------------*/
void srmig3_close(ssr3d ssr,
		  tap3d tap,
		  slo3d s_s,
		  slo3d s_r,
		  sroperator3d srop)
/*< free allocated storage >*/
{
    ssr3_close(ssr,tap);
    
    slow3_close(s_s);
    slow3_close(s_r);

    free(**srop->ww_s); free( *srop->ww_s); free( srop->ww_s);
    free(**srop->ww_r); free( *srop->ww_r); free( srop->ww_r);
}

/*------------------------------------------------------------*/
void srmig3(fslice sdat /* source   data [nw][ny][nx] */,
	    fslice rdat /* receiver data [nw][ny][nx] */,
	    fslice imag /*         image [nz][ny][nx] */,
	    fslice cigs,
	    void (*imop)(int,int),
	    int ompchunk,
	    sroperator3d srop
    )
/*< Apply S/R migration >*/
{
    int imz,iw,ie;
    sf_complex ws,wr;
    
    int ompith=0;

    for (ie=0; ie<srop->ae.n; ie++) {

#ifdef _OPENMP
#pragma omp parallel for schedule(static)   \
    private(ompith,iw,ws,wr,imz)	    \
    shared(sdat,rdat,ie,srop)
#endif
	for (iw=0; iw<srop->aw.n; iw++) {
	    ompith=omp_get_thread_num();
#ifdef _OPENMP	    
#pragma omp critical
	    if(srop->verb) sf_warning ("(ith=%d) ... <iw=%3d of %3d> ... <ie=%3d of %3d>",
				 ompith,iw+1,srop->aw.n,ie+1,srop->ae.n);
#endif
	    
	    ws = sf_cmplx(srop->eps*srop->aw.d,+(srop->aw.o+iw*srop->aw.d)); //      causal
	    wr = sf_cmplx(srop->eps*srop->aw.d,-(srop->aw.o+iw*srop->aw.d)); // anti-causal
	    
#ifdef _OPENMP	    
#pragma omp critical
	    {
		fslice_get(sdat,ie*srop->aw.n+iw,srop->ww_s[ompith][0]);
		fslice_get(rdat,ie*srop->aw.n+iw,srop->ww_r[ompith][0]);
	    }
#endif	    
	    taper2d(srop->ww_s[ompith],srop->tap);
	    taper2d(srop->ww_r[ompith],srop->tap);	    
	    
	    fslice_get(srop->s_s->slice, 0, srop->s_s->so[ompith][0]);
	    fslice_get(srop->s_r->slice, 0, srop->s_r->so[ompith][0]);

	    for (imz=0; imz<srop->amz.n-1; imz++) {

		fslice_get(srop->s_s->slice, imz+1, srop->s_s->ss[ompith][0]);
		fslice_get(srop->s_r->slice, imz+1, srop->s_r->ss[ompith][0]);

		ssr3_ssf(ws,srop->ww_s[ompith],srop->ssr,srop->tap,srop->s_s,imz,ompith);
		ssr3_ssf(wr,srop->ww_r[ompith],srop->ssr,srop->tap,srop->s_r,imz,ompith);
		
		slow3_advance(srop->s_s,ompith);
		slow3_advance(srop->s_r,ompith);

		img3store(imz,srop->ww_s,srop->ww_r,ompith);
		
	    } // z 

	    imop(iw,ompith); // imaging operator
	    
	} // w
	
    } // e
}

