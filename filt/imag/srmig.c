/* 3-D SSR migration/modeling using extended split-step */

/*
  Copyright (C) 2006 Colorado School of Mines
  Copyright (C) 2004 University of Texas at Austin
  
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

#include "srmig.h"
#include "taper.h"
#include "slowref.h"
#include "ssr.h"
#include "img.h"

#include "slice.h"
/*^*/

#define SOOP(a) for(ily=0;ily<aly.n;ily++){ \
                for(ilx=0;ilx<alx.n;ilx++){ {a} }}

static sf_axa aw,alx,aly,amx,amy,amz,ae;
static bool verb;
static float eps;

static sf_complex **ww_s,**ww_r;

static int            *nr_s, *nr_r, *nr; /* number of references */
static float         **sm_s,**sm_r,**sm; /* ref slo squared */
static fslice        slow_s,slow_r,slow; /* slowness slice */
static float         **ss_s,**ss_r,**ss; /* slowness */
static float         **so_s,**so_r,**so; /* slowness */

void srmig_init(bool verb_,
		float eps_,
		float dtmax,
		sf_axis ae_      /* experiments (e.g. shots) */,
		sf_axis aw_      /* frequency */,
		sf_axis amx_     /* i-line (data) */,
		sf_axis amy_     /* x-line (data) */,
		sf_axis amz_     /* depth */,
		sf_axis alx_     /* i-line (slowness/image) */,
		sf_axis aly_     /* x-line (slowness/image) */,
		int tmx, int tmy /* taper size */,
		int pmx, int pmy /* padding in the k domain */
    )
/*< initialize SR migration >*/
{
    float dsmax;

    verb=verb_;
    eps = eps_;

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
    ssr_init(amz_ ,
	     amx_,amy_,
	     alx_,aly_,
	     pmx ,pmy,
	     tmx ,tmy,
	     dsmax);

    /* precompute taper */
    taper2_init(amy.n,
		amx.n,
		SF_MIN(tmy,amy.n-1),
		SF_MIN(tmx,amx.n-1),
		true,
		true);

    /* allocate wavefield storage */
    ww_s = sf_complexalloc2(amx.n,amy.n);
    ww_r = sf_complexalloc2(amx.n,amy.n);
}
/*------------------------------------------------------------*/

void srmig_pw_init(float  dtmax,
		   int    nrmax,   /* maximum number of references */
		   fslice slow_)
/*< initialize P-wave slowness >*/
{
    int   imz, jj;
    float dsmax;

    dsmax = dtmax/amz.d;

    /*------------------------------------------------------------*/
    /* slowness: downgoing wavefield */
    ss = sf_floatalloc2(alx.n,aly.n);  /* slowness */
    so = sf_floatalloc2(alx.n,aly.n);  /* slowness */
    sm = sf_floatalloc2(nrmax,amz.n); /* ref slowness squared */
    nr = sf_intalloc         (amz.n); /* nr of ref slownesses */
    slow = slow_;
    for (imz=0; imz<amz.n; imz++) {
	fslice_get(slow,imz,ss[0]);
	
	nr[imz] = slowref(nrmax,dsmax,alx.n*aly.n,ss[0],sm[imz]);
	if (verb) sf_warning("nr[%d]=%d",imz,nr[imz]);
    }
    for (imz=0; imz<amz.n-1; imz++) {
	for (jj=0; jj<nr[imz]; jj++) {
	    sm[imz][jj] = 0.5*(sm[imz][jj]+sm[imz+1][jj]);
	}
    }
}

/*------------------------------------------------------------*/

void srmig_cw_init(float  dtmax,
		   int    nrmax,   /* maximum number of references */
		   fslice slow_s_,
		   fslice slow_r_)
/*< initialize C-wave slowness >*/
{
    int   imz, jj;
    float dsmax;

    dsmax = dtmax/amz.d;

    /*------------------------------------------------------------*/
    /* slowness: downgoing wavefield */
    ss_s = sf_floatalloc2(alx.n,aly.n);  /* slowness */
    so_s = sf_floatalloc2(alx.n,aly.n);  /* slowness */
    sm_s = sf_floatalloc2(nrmax,amz.n); /* ref slowness squared */
    nr_s = sf_intalloc         (amz.n); /* nr of ref slownesses */
    slow_s = slow_s_;
    for (imz=0; imz<amz.n; imz++) {
	fslice_get(slow_s,imz,ss_s[0]);
	
	nr_s[imz] = slowref(nrmax,dsmax,alx.n*aly.n,ss_s[0],sm_s[imz]);
	if (verb) sf_warning("nr_s[%d]=%d",imz,nr_s[imz]);
    }
    for (imz=0; imz<amz.n-1; imz++) {
	for (jj=0; jj<nr_s[imz]; jj++) {
	    sm_s[imz][jj] = 0.5*(sm_s[imz][jj]+sm_s[imz+1][jj]);
	}
    }
    /*------------------------------------------------------------*/
    /* slowness:   upgoing wavefield */
    ss_r = sf_floatalloc2(alx.n,aly.n);  /* slowness */
    so_r = sf_floatalloc2(alx.n,aly.n);  /* slowness */
    sm_r = sf_floatalloc2(nrmax,amz.n); /* ref slowness squared */
    nr_r = sf_intalloc         (amz.n); /* nr of ref slownesses */
    slow_r = slow_r_;
    for (imz=0; imz<amz.n; imz++) {
	fslice_get(slow_r,imz,ss_r[0]);
	
	nr_r[imz] = slowref(nrmax,dsmax,alx.n*aly.n,ss_r[0],sm_r[imz]);
	if (verb) sf_warning("nr_r[%d]=%d",imz,nr_r[imz]);
    }
    for (imz=0; imz<amz.n-1; imz++) {
	for (jj=0; jj<nr_r[imz]; jj++) {
	    sm_r[imz][jj] = 0.5*(sm_r[imz][jj]+sm_r[imz+1][jj]);
	}
    }
    /*------------------------------------------------------------*/
}

/*------------------------------------------------------------*/

void srmig_pw_close(void)
/*< free slowness storage (P waves) >*/
{
    free( *ss); free( ss);
    free( *so); free( so);
    free( *sm); free( sm);
    ;           free( nr);
}

/*------------------------------------------------------------*/

void srmig_cw_close(void)
/*< free slowness storage (C waves) >*/
{
    free( *ss_s); free( ss_s);
    free( *so_s); free( so_s);
    free( *sm_s); free( sm_s);
    ;             free( nr_s);

    free( *ss_r); free( ss_r);
    free( *so_r); free( so_r);
    free( *sm_r); free( sm_r);
    ;             free( nr_r);
}

/*------------------------------------------------------------*/

void srmig_close(void)
/*< free allocated storage >*/
{
    ssr_close();
    
    free( *ww_s); free( ww_s);
    free( *ww_r); free( ww_r);
}

/*------------------------------------------------------------*/
void srmig_pw(fslice sdat /* source   data [nw][ny][nx] */,
	      fslice rdat /* receiver data [nw][ny][nx] */,
	      fslice imag /*         image [nz][ny][nx] */,
	      void (*imop)( fslice, int)
    )
/*< Apply S/R migration >*/
{
    int imz,iw,ilx,ily,ie;
    sf_complex ws,wr;
    
    for (ie=0; ie<ae.n; ie++) {
	for (iw=0; iw<aw.n; iw++) {
	    if(verb) sf_warning ("iw=%3d of %3d:   ie=%3d of %3d",
				 iw+1,aw.n,ie+1,ae.n);
	    
	    ws = sf_cmplx(eps*aw.d,+(aw.o+iw*aw.d)); /*      causal */
	    wr = sf_cmplx(eps*aw.d,-(aw.o+iw*aw.d)); /* anti-causal */
	    
	    fslice_get(sdat,ie*aw.n+iw,ww_s[0]); taper2(ww_s);
	    fslice_get(rdat,ie*aw.n+iw,ww_r[0]); taper2(ww_r);

	    fslice_get(slow,0,so[0]);
	    for (imz=0; imz<amz.n-1; imz++) {
		fslice_get(slow,imz+1,ss[0]);

		ssr_ssf(ws,ww_s,so,ss,nr[imz],sm[imz]);
		ssr_ssf(wr,ww_r,so,ss,nr[imz],sm[imz]);
		
		SOOP( so[ily][ilx] = ss[ily][ilx]; );
		
		imgstore(imz,ww_s,ww_r);
	    } /* z */
	    
	    imop(imag,iw);
	} /* w */
    } /* e */
}

/*------------------------------------------------------------*/

void srmig_cw(fslice sdat /* source   data [nw][ny][nx] */,
	      fslice rdat /* receiver data [nw][ny][nx] */,
	      fslice imag /*         image [nz][ny][nx] */,
	      void (*imop)( fslice,int)
    )
/*< Apply S/R migration >*/
{
    int imz,iw,ilx,ily,ie;
    sf_complex ws,wr;
    
    for (ie=0; ie<ae.n; ie++) {
	for (iw=0; iw<aw.n; iw++) {
	    if(verb) sf_warning ("iw=%3d of %3d:   ie=%3d of %3d",
				 iw+1,aw.n,ie+1,ae.n);
	    
	    ws = sf_cmplx(eps*aw.d,+(aw.o+iw*aw.d)); /*      causal */
	    wr = sf_cmplx(eps*aw.d,-(aw.o+iw*aw.d)); /* anti-causal */
	    
	    fslice_get(sdat,ie*aw.n+iw,ww_s[0]); taper2(ww_s);
	    fslice_get(rdat,ie*aw.n+iw,ww_r[0]); taper2(ww_r);
	    
	    fslice_get(slow_s,0,so_s[0]);
	    fslice_get(slow_r,0,so_r[0]);
	    for (imz=0; imz<amz.n-1; imz++) {
		fslice_get(slow_s,imz+1,ss_s[0]);
		fslice_get(slow_r,imz+1,ss_r[0]);

		ssr_ssf(ws,ww_s,so_s,ss_s,nr_s[imz],sm_s[imz]);
		ssr_ssf(wr,ww_r,so_r,ss_r,nr_r[imz],sm_r[imz]);
		
		SOOP( so_s[ily][ilx] = ss_s[ily][ilx]; );
		SOOP( so_r[ily][ilx] = ss_r[ily][ilx]; );

		imgstore(imz,ww_s,ww_r);
	    } /* z */

	    imop(imag,iw);
	} /* w */
    } /* e */
}

/*------------------------------------------------------------*/

