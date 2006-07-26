/* 3-D shot-record modeling using extended split-step */
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

#include "srmod.h"
#include "taper.h"
#include "slowref.h"
#include "ssr.h"

#include "slice.h"
/*^*/

#define LOOP(a) for( iy=0; iy< ay.n; iy++){ \
                for( ix=0; ix< ax.n; ix++){ {a} }}
#define SOOP(a) for(ily=0;ily<aly.n;ily++){ \
                for(ilx=0;ilx<alx.n;ilx++){ {a} }}

static sf_axa az,aw,alx,aly,ax,ay;
static bool verb, incore;
static float eps;

static float         **ss; /* slowness */
static float         **so; /* slowness */

static int            *nr_s, *nr_r, *nr; /* number of refs */
static float         **sm_s,**sm_r,**sm; /* ref slo squared */
static fslice        slow_s,slow_r,slow; /* slowness slice */
static fslice        wtmp;

/* read one depth level in memory */
static float         **rr;
static sf_complex **ww_s;
static sf_complex **ww_r;

/* read all depth levels in memory */
static float         ***rrr;
static sf_complex ***www;
static float         ***ssp;
static float         ***sss;

void srmod_init(bool verb_,
		bool incore_,  /* keep shot wavefield in core */
		float eps_,
		float dtmax,
		sf_axis az_    /* depth */,
		sf_axis aw_    /* frequency */,
		sf_axis ax_    /* i-line (data) */,
		sf_axis ay_    /* x-line (data) */,
		sf_axis alx_   /* i-line (slowness/image) */,
		sf_axis aly_   /* x-line (slowness/image) */,
		int tx, int ty /* taper size */,
		int px, int py /* padding in the k domain */
    )
/*< initialize S/R modeling >*/
{
    float dsmax;

    verb   = verb_;
    incore = incore_;
    eps    = eps_;
    
    az = sf_nod(az_);
    aw = sf_nod(aw_);
    ax = sf_nod(ax_);
    ay = sf_nod(ay_);
    alx= sf_nod(alx_);
    aly= sf_nod(aly_);

    /* from hertz to radian */
    aw.d *= 2.*SF_PI; 
    aw.o *= 2.*SF_PI;

    dsmax  = dtmax/az.d;

    /* SSR */
    ssr_init(az_ ,
	     ax_ ,ay_,
	     alx_,aly_,
	     px  ,py,
	     tx  ,ty,
	     dsmax);

    /* precompute taper */
    taper2_init(ay.n,
		ax.n,
		SF_MIN(ty,ay.n-1),
		SF_MIN(tx,ax.n-1),
		true,
		true);

    /* compute reference slowness */
    ss = sf_floatalloc2(alx.n,aly.n); /* slowness */
    so = sf_floatalloc2(alx.n,aly.n); /* slowness */

    ww_s = sf_complexalloc2(ax.n,ay.n); /* wavefield storage */
    ww_r = sf_complexalloc2(ax.n,ay.n);
    rr   = sf_floatalloc2  (ax.n,ay.n); /* reflectivity storage */
    
    if(incore) {
	www = sf_complexalloc3(ax.n,ay.n,az.n);
	rrr = sf_floatalloc3  (ax.n,ay.n,az.n);
    } else {
	wtmp  = fslice_init( ax.n * ay.n, az.n,sizeof(sf_complex));
    }
}

void srmod_pw_init( float dtmax,
		    int   nrmax      /* maximum number of references */,
		    fslice slow_)
/*< initialize P-wave S/R modeling >*/
{
    int   iz, jj;
    float dsmax;

    dsmax  = dtmax/az.d;

    sm= sf_floatalloc2 (nrmax,az.n); /* ref slowness squared*/
    nr= sf_intalloc          (az.n); /* nr of ref slownesses */
    slow= slow_;
    for (iz=0; iz<az.n; iz++) {
	fslice_get(slow,iz,ss[0]);
	
	nr[iz] = slowref(nrmax,dsmax,alx.n*aly.n,ss[0],sm[iz]);
	if (verb) sf_warning("nr[%d]=%d",iz,nr[iz]);
    }
    for (iz=0; iz<az.n-1; iz++) {
	for (jj=0; jj<nr[iz]; jj++) {
	    sm[iz][jj] = 0.5*(sm[iz][jj]+sm[iz+1][jj]);
	}
    }

    if(incore) {
	ssp = sf_floatalloc3  (alx.n,aly.n,az.n);
    }

}

void srmod_cw_init( float dtmax,
		    int   nrmax      /* maximum number of references */,
		    fslice slow_s_,
		    fslice slow_r_)
/*< initialize C-wave S/R modeling >*/
{
    int   iz, jj;
    float dsmax;

    dsmax  = dtmax/az.d;

    /*------------------------------------------------------------*/
    /* slowness: downgoing wavefield */
    sm_s= sf_floatalloc2 (nrmax,az.n); /* ref slowness squared*/
    nr_s= sf_intalloc          (az.n); /* nr of ref slownesses */
    slow_s= slow_s_;
    for (iz=0; iz<az.n; iz++) {
	fslice_get(slow_s,iz,ss[0]);
	
	nr_s[iz] = slowref(nrmax,dsmax,alx.n*aly.n,ss[0],sm_s[iz]);
	if (verb) sf_warning("nr_s[%d]=%d",iz,nr_s[iz]);
    }
    for (iz=0; iz<az.n-1; iz++) {
	for (jj=0; jj<nr_s[iz]; jj++) {
	    sm_s[iz][jj] = 0.5*(sm_s[iz][jj]+sm_s[iz+1][jj]);
	}
    }
    /*------------------------------------------------------------*/
    /* slowness: up-going wavefield */
    sm_r= sf_floatalloc2 (nrmax,az.n); /* ref slowness squared*/
    nr_r= sf_intalloc          (az.n); /* nr of ref slownesses */
    slow_r= slow_r_;
    for (iz=0; iz<az.n; iz++) {
	fslice_get(slow_r,iz,ss[0]);
	
	nr_r[iz] = slowref(nrmax,dsmax,alx.n*aly.n,ss[0],sm_r[iz]);
	if (verb) sf_warning("nr_r[%d]=%d",iz,nr_r[iz]);
    }
    for (iz=0; iz<az.n-1; iz++) {
	for (jj=0; jj<nr_r[iz]; jj++) {
	    sm_r[iz][jj] = 0.5*(sm_r[iz][jj]+sm_r[iz+1][jj]);
	}
    }
    /*------------------------------------------------------------*/

    if(incore) {
	ssp = sf_floatalloc3  (alx.n,aly.n,az.n);
	sss = sf_floatalloc3  (alx.n,aly.n,az.n);
    }
}


/*------------------------------------------------------------*/
void srmod_pw_close(void)
/*< free P-wave slowness storage >*/
{
    free( *sm); free( sm);
    ;           free( nr);

    if(incore) {
	free(**ssp); free(*ssp); free(ssp);
    }
}

void srmod_cw_close(void)
/*< free C-wave slowness storage >*/
{
    free( *sm_s); free( sm_s);
    ;             free( nr_s);
    free( *sm_r); free( sm_r);
    ;             free( nr_r);

    if(incore) {
	free(**ssp); free(*ssp); free(ssp);
	free(**sss); free(*sss); free(sss);
    }
}

void srmod_close(void)
/*< free allocated storage >*/
{
    ssr_close();

    free( *ss); free( ss);
    free( *so); free( so);
    
    free( *ww_s); free( ww_s);
    free( *ww_r); free( ww_r);
    free( *rr  ); free( rr  );

    if(incore) {
	free(**www); free( *www); free( www);
	free(**rrr); free( *rrr); free( rrr);
    } else {
	fslice_close(wtmp);
    }
}

/*------------------------------------------------------------*/
void srmod_pw(fslice dwfl /* source   data [nw][ny][nx] */,
	      fslice uwfl /* receiver data [nw][ny][nx] */,
	      fslice refl
    )
/*< Apply P-wave S/R modeling >*/
{
    int iz,iw,ix,iy,ilx,ily;
    sf_complex w;

    if(incore) {

	/* read reflectivity */
	for( iz=0; iz<az.n; iz++) {
	    fslice_get(refl,iz,rr[0]);
	    LOOP( rrr[iz][ iy][ ix] = rr[iy][ix];);

	    fslice_get(slow,iz,so[0]);
	    SOOP( ssp[iz][ily][ilx] = so[ily][ilx];);
	}

	for (iw=0; iw<aw.n; iw++) {
	    if(verb) sf_warning("iw=%3d of %3d (in core)",iw+1,aw.n);
	    w = sf_cmplx(eps*aw.d,aw.o+iw*aw.d);
	    
	    /* downgoing wavefield */
	    fslice_get(dwfl,iw,ww_s[0]); taper2(ww_s);
	    LOOP( www[0][iy][ix] = ww_s[iy][ix]; );
	    for (iz=0; iz<az.n-1; iz++) {
		ssr_ssf(w,ww_s,ssp[iz],ssp[iz+1],nr[iz],sm[iz]);
		LOOP( www[iz+1][iy][ix] = ww_s[iy][ix]; );
	    }
	    
	    /* upgoing wavefield */
	    LOOP( ww_r[iy][ix] = sf_cmplx(0.,0.); );	    
	    for (iz=az.n-1; iz>0; iz--) {
#ifdef SF_HAS_COMPLEX_H
		LOOP( ww_r[iy][ix] += www[iz][iy][ix]*rrr[iz][iy][ix]; );
#else
		LOOP( ww_r[iy][ix] = sf_cadd(ww_r[iy][ix],
					     sf_crmul(www[iz][iy][ix],
						      rrr[iz][iy][ix])); );
#endif
		ssr_ssf(w,ww_r,ssp[iz],ssp[iz-1],nr[iz-1],sm[iz-1]);
	    }

	    fslice_put(uwfl,iw,ww_r[0]);
	} /* iw */

    } else {

	for (iw=0; iw<aw.n; iw++) {
	    if(verb) sf_warning("iw=%3d of %3d",iw+1,aw.n);
	    w = sf_cmplx(eps*aw.d,aw.o+iw*aw.d);
	    
	    /* downgoing wavefield */
	    fslice_get(dwfl,iw,ww_s[0]); taper2(ww_s);
	    fslice_put(wtmp, 0,ww_s[0]);
	    
	    fslice_get(slow,0,so[0]);
	    for (iz=0; iz<az.n-1; iz++) {
		fslice_get(slow,iz+1,ss[0]);
		ssr_ssf(w,ww_s,so,ss,nr[iz],sm[iz]);
		SOOP( so[ily][ilx] = ss[ily][ilx]; );
		
		fslice_put(wtmp,iz+1,ww_s[0]);
	    }
	    
	    /* upgoing wavefield */
	    LOOP( ww_r[iy][ix] = sf_cmplx(0.,0.); );
	    
	    fslice_get(slow,az.n-1,so[0]);
	    for (iz=az.n-1; iz>0; iz--) {
		fslice_get(slow,iz-1,ss[0]);
		
		fslice_get(wtmp,iz,ww_s[0]); 
		fslice_get(refl,iz,rr[0]); /* reflectivity */

#ifdef SF_HAS_COMPLEX_H		
		LOOP( ww_s[iy][ix] *= rr  [iy][ix];
		      ww_r[iy][ix] += ww_s[iy][ix]; );
#else
		LOOP( ww_s[iy][ix] = sf_crmul(ww_s[iy][ix],rr[iy][ix]);
		      ww_r[iy][ix] = sf_cadd(ww_r[iy][ix],ww_s[iy][ix]); );
#endif		
		ssr_ssf(w,ww_r,so,ss,nr[iz-1],sm[iz-1]);
		SOOP( so[ily][ilx] = ss[ily][ilx]; );
	    }
	    fslice_put(uwfl,iw,ww_r[0]);
	} /* iw */
    } /* in core */
}

/*------------------------------------------------------------*/

void srmod_cw(fslice dwfl /* source   data [nw][ny][nx] */,
	      fslice uwfl /* receiver data [nw][ny][nx] */,
	      fslice refl
    )
/*< Apply C-wave S/R modeling >*/
{
    int iz,iw,ix,iy,ilx,ily;
    sf_complex w;

    if(incore) {

	/* read reflectivity */
	for( iz=0; iz<az.n; iz++) {
	    fslice_get(refl,iz,rr[0]);
	    LOOP( rrr[iz][iy][ix] = rr[iy][ix];);

	    fslice_get(slow_s,iz,so[0]);
	    SOOP( ssp[iz][ily][ilx] = so[ily][ilx];);

	    fslice_get(slow_r,iz,so[0]);
	    SOOP( sss[iz][ily][ilx] = so[ily][ilx];);
	}
	
	for (iw=0; iw<aw.n; iw++) {
	    if(verb) sf_warning("iw=%3d of %3d (in core)",iw+1,aw.n);
	    w = sf_cmplx(eps*aw.d,aw.o+iw*aw.d);
	    
	    /* downgoing wavefield */
	    fslice_get(dwfl,iw,ww_s[0]); taper2(ww_s);
	    LOOP( www[0][iy][ix] = ww_s[iy][ix]; );	    
	    for (iz=0; iz<az.n-1; iz++) {
		ssr_ssf(w,ww_s,ssp[iz],ssp[iz+1],nr_s[iz],sm_s[iz]);
		LOOP( www[iz+1][iy][ix] = ww_s[iy][ix]; );
	    }
	    
	    /* upgoing wavefield */
	    LOOP( ww_r[iy][ix] = sf_cmplx(0.,0.); );
	    for (iz=az.n-1; iz>0; iz--) {
#ifdef SF_HAS_COMPLEX_H
		LOOP( ww_r[iy][ix] += www[iz][iy][ix]*rrr[iz][iy][ix]; );
#else
		LOOP( ww_r[iy][ix] = sf_cadd(ww_r[iy][ix],
					     sf_crmul(www[iz][iy][ix],
						      rrr[iz][iy][ix])); );
#endif
		ssr_ssf(w,ww_r,sss[iz],sss[iz-1],nr[iz],sm[iz]);
	    }
	    fslice_put(uwfl,iw,ww_r[0]);
	} /* iw */

    } else {

	for (iw=0; iw<aw.n; iw++) {
	    if(verb) sf_warning("iw=%3d of %3d",iw+1,aw.n);
	    w = sf_cmplx(eps*aw.d,aw.o+iw*aw.d);
	    
	    /* downgoing wavefield */
	    fslice_get(dwfl,iw,ww_s[0]); taper2(ww_s);
	    
	    fslice_put(wtmp,0,ww_s[0]);
	    
	    fslice_get(slow_s,0,so[0]);
	    for (iz=0; iz<az.n-1; iz++) {
		fslice_get(slow_s,iz+1,ss[0]);
		ssr_ssf(w,ww_s,so,ss,nr_s[iz],sm_s[iz]);
		SOOP( so[ily][ilx] = ss[ily][ilx]; );
		
		fslice_put(wtmp,iz+1,ww_s[0]);
	    }
	    
	    /* upgoing wavefield */
	    LOOP( ww_r[iy][ix] = sf_cmplx(0.,0.); );
	    
	    fslice_get(slow_r,az.n-1,so[0]);
	    for (iz=az.n-1; iz>0; iz--) {
		fslice_get(slow_r,iz-1,ss[0]);
		
		fslice_get(wtmp,iz,ww_s[0]); 
		fslice_get(refl,iz,rr[0]); /* reflectivity */
#ifdef SF_HAS_COMPLEX_H
		LOOP( ww_s[iy][ix] *= rr  [iy][ix];
		      ww_r[iy][ix] += ww_s[iy][ix]; );
#else
		LOOP( ww_s[iy][ix] = sf_crmul(ww_s[iy][ix],rr[iy][ix]);
		      ww_r[iy][ix] = sf_cadd(ww_r[iy][ix],ww_s[iy][ix]); );
#endif
		
		ssr_ssf(w,ww_r,so,ss,nr_r[iz],sm_r[iz]);
		SOOP( so[ily][ilx] = ss[ily][ilx]; );
	    }
	    fslice_put(uwfl,iw,ww_r[0]);
	} /* iw */
    } /* in core */

}
