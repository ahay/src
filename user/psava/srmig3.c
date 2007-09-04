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
#include "slowref.h"
#include "ssr3.h"
#include "img3.h"

#include "slice.h"
/*^*/

#include "weutil.h"
/*^*/

#define SOOP(a) for(ily=0;ily<aly.n;ily++){ \
                for(ilx=0;ilx<alx.n;ilx++){ {a} }}

static sf_axa aw,alx,aly,amx,amy,amz,ae;
static bool verb;
static float eps;

static fslice     slow_s, slow_r; /* slowness slice */
static sf_complex***ww_s,***ww_r;
static int         *nr_s,  *nr_r; /* number of references */
static float      **sm_s, **sm_r; /* ref slo squared */

static float     ***ss_s,***ss_r; /* slowness */
static float     ***so_s,***so_r; /* slowness */

static int ompnth;

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
    int   imz, jj;

    sroperator3d srop;
    srop = (sroperator3d) sf_alloc(1,sizeof(*srop));
    
    verb=verb_;
    eps = eps_;

    ompnth = ompnth_;

    aw  = sf_nod(aw_);
    amx = sf_nod(amx_);
    amy = sf_nod(amy_);
    amz = sf_nod(amz_);
    alx = sf_nod(alx_);
    aly = sf_nod(aly_);
    ae  = sf_nod(ae_);

    /* from hertz to radian */
    aw.d *= 2.*SF_PI; 
    aw.o *= 2.*SF_PI;

    dsmax = dtmax/amz.d;

    /* SSR */
    srop->ssr = ssr3_init(amz_ ,
			  amx_,amy_,
			  alx_,aly_,
			  pmx ,pmy,
			  tmx ,tmy,
			  dsmax,
			  ompnth);
    
    /* taper */
    srop->tap = taper_init(amx.n,amy.n,1,
			   SF_MIN(tmx,amx.n-1),SF_MIN(tmy,amy.n-1),0,
			   true,true,false);
    
    /* allocate wavefield storage */
    ww_s = sf_complexalloc3(amx.n,amy.n,ompnth);
    ww_r = sf_complexalloc3(amx.n,amy.n,ompnth);
    
    /*------------------------------------------------------------*/
    /* slowness: downgoing wavefield */
    ss_s = sf_floatalloc3(alx.n,aly.n,ompnth);  /* slowness */
    so_s = sf_floatalloc3(alx.n,aly.n,ompnth);  /* slowness */
    sm_s = sf_floatalloc2(nrmax,amz.n);  /* ref slowness squared */
    nr_s = sf_intalloc         (amz.n);  /* nr of ref slownesses */
    slow_s = slow_s_;
    for (imz=0; imz<amz.n; imz++) {
	fslice_get(slow_s,imz,ss_s[0][0]);
	
	nr_s[imz] = slowref(nrmax,dsmax,alx.n*aly.n,ss_s[0][0],sm_s[imz]);
    }
    for (imz=0; imz<amz.n-1; imz++) {
	for (jj=0; jj<nr_s[imz]; jj++) {
	    sm_s[imz][jj] = 0.5*(sm_s[imz][jj]+sm_s[imz+1][jj]);
	}
    }
    /*------------------------------------------------------------*/
    /* slowness:   upgoing wavefield */
    ss_r = sf_floatalloc3(alx.n,aly.n,ompnth);  /* slowness */
    so_r = sf_floatalloc3(alx.n,aly.n,ompnth);  /* slowness */
    sm_r = sf_floatalloc2(nrmax,amz.n);  /* ref slowness squared */
    nr_r = sf_intalloc         (amz.n);  /* nr of ref slownesses */
    slow_r = slow_r_;
    for (imz=0; imz<amz.n; imz++) {
	fslice_get(slow_r,imz,ss_r[0][0]);
	
	nr_r[imz] = slowref(nrmax,dsmax,alx.n*aly.n,ss_r[0][0],sm_r[imz]);
    }
    for (imz=0; imz<amz.n-1; imz++) {
	for (jj=0; jj<nr_r[imz]; jj++) {
	    sm_r[imz][jj] = 0.5*(sm_r[imz][jj]+sm_r[imz+1][jj]);
	}
    }

    return srop;
}

/*------------------------------------------------------------*/
void srmig3_close(ssr3d ssr,
		  tap3d tap)
/*< free allocated storage >*/
{
    ssr3_close(ssr,tap);
    
    free(**ss_s); free( *ss_s); free( ss_s);
    free(**so_s); free( *so_s); free( so_s);
    free( *sm_s); free( sm_s);
    ;             free( nr_s);

    free(**ss_r); free( *ss_r); free( ss_r);
    free(**so_r); free( *so_r); free( so_r);
    free( *sm_r); free( sm_r);
    ;             free( nr_r);

    free(**ww_s); free( *ww_s); free( ww_s);
    free(**ww_r); free( *ww_r); free( ww_r);
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
    int imz,iw,ilx,ily,ie;
    sf_complex ws,wr;
    
    int ompith=0;

    for (ie=0; ie<ae.n; ie++) {

#ifdef _OPENMP
#pragma omp parallel for schedule(static)				\
    private(ompith,iw,ws,wr,imz,ilx,ily)				\
    shared(aw,eps,sdat,rdat,ie,ww_s,ww_r,slow_s,slow_r,so_s,so_r,ss_s,ss_r,amz,nr_s,nr_r,sm_s,sm_r,alx,aly)
#endif
	for (iw=0; iw<aw.n; iw++) {
	    ompith=omp_get_thread_num();
#ifdef _OPENMP	    
#pragma omp critical
	    if(verb) sf_warning ("(ith=%d) ... <iw=%3d of %3d> ... <ie=%3d of %3d>",
				 ompith,iw+1,aw.n,ie+1,ae.n);
#endif
	    
	    ws = sf_cmplx(eps*aw.d,+(aw.o+iw*aw.d)); //      causal
	    wr = sf_cmplx(eps*aw.d,-(aw.o+iw*aw.d)); // anti-causal
	    
#ifdef _OPENMP	    
#pragma omp critical
	    {
		fslice_get(sdat,ie*aw.n+iw,ww_s[ompith][0]);
		fslice_get(rdat,ie*aw.n+iw,ww_r[ompith][0]);
	    }
#endif	    
	    taper2d(ww_s[ompith],srop->tap);
	    taper2d(ww_r[ompith],srop->tap);	    
	    
	    fslice_get(slow_s,0,so_s[ompith][0]);
	    fslice_get(slow_r,0,so_r[ompith][0]);
	    for (imz=0; imz<amz.n-1; imz++) {
		fslice_get(slow_s,imz+1,ss_s[ompith][0]);
		fslice_get(slow_r,imz+1,ss_r[ompith][0]);

		ssr3_ssf(ws,ww_s[ompith],so_s[ompith],ss_s[ompith],nr_s[imz],sm_s[imz],ompith,srop->ssr,srop->tap);
		ssr3_ssf(wr,ww_r[ompith],so_r[ompith],ss_r[ompith],nr_r[imz],sm_r[imz],ompith,srop->ssr,srop->tap);
		
		SOOP( so_s[ompith][ily][ilx] = ss_s[ompith][ily][ilx]; );
		SOOP( so_r[ompith][ily][ilx] = ss_r[ompith][ily][ilx]; );
		
		img3store(imz,ww_s,ww_r,ompith);
		
	    } // z 
	    imop(iw,ompith);
	    
	} // w
	
    } // e
}

