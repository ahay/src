/* 3-D shot-record modeling using extended split-step */
/*
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

static axa az,aw,alx,aly,ax,ay,ae;
static bool verb;
static float eps;

static float         **qq;
static float complex **ud;
static float complex **uu;
static float         **rr; /* reflectivity */

static float         **ss; /* slowness */
static float         **so; /* slowness */

static int            *nrd, *nru; /* number of refs */
static float         **smd,**smu; /* ref slo squared */
static fslice        slowd,slowu; /* slowness slice */

void srmod_init(bool verb_,
		float eps_,
		float dtmax,
		axa ae_        /* experiments (e.g. shots) */,
		axa az_        /* depth */,
		axa aw_        /* frequency */,
		axa ax_        /* i-line (data) */,
		axa ay_        /* x-line (data) */,
		axa alx_       /* i-line (slowness/image) */,
		axa aly_       /* x-line (slowness/image) */,
		int tx, int ty /* taper size */,
		int px, int py /* padding in the k domain */,
		int nrmax      /* maximum number of references */,
		fslice slowd_,
		fslice slowu_)
/*< initialize >*/
{
    int   iz, jj;
    float ds;

    verb=verb_;
    eps = eps_;
    
    az = az_;
    aw = aw_;
    ax = ax_;
    ay = ay_;
    alx= alx_;
    aly= aly_;
    ae = ae_;

    /* from hertz to radian */
    aw.d *= 2.*SF_PI; 
    aw.o *= 2.*SF_PI;

    ds  = dtmax/az.d;

    /* SSR */
    ssr_init(az_ ,
	     ax_ ,ay_,
	     alx_,aly_,
	     px  ,py,
	     tx  ,ty,
	     ds);

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


    /*------------------------------------------------------------*/
    /* slowness: downgoing wavefield */
    smd= sf_floatalloc2 (nrmax,az.n); /* ref slowness squared*/
    nrd= sf_intalloc          (az.n); /* nr of ref slownesses */
    slowd= slowd_;
    for (iz=0; iz<az.n; iz++) {
	fslice_get(slowd,iz,ss[0]);
	
	nrd[iz] = slowref(nrmax,ds,alx.n*aly.n,ss[0],smd[iz]);
	if (verb) sf_warning("nrd[%d]=%d",iz,nrd[iz]);
    }
    for (iz=0; iz<az.n-1; iz++) {
	for (jj=0; jj<nrd[iz]; jj++) {
	    smd[iz][jj] = 0.5*(smd[iz][jj]+smd[iz+1][jj]);
	}
    }
    /*------------------------------------------------------------*/
    /* slowness: up-going wavefield */
    smu= sf_floatalloc2 (nrmax,az.n); /* ref slowness squared*/
    nru= sf_intalloc          (az.n); /* nr of ref slownesses */
    slowu= slowu_;
    for (iz=0; iz<az.n; iz++) {
	fslice_get(slowu,iz,ss[0]);
	
	nru[iz] = slowref(nrmax,ds,alx.n*aly.n,ss[0],smu[iz]);
	if (verb) sf_warning("nru[%d]=%d",iz,nru[iz]);
    }
    for (iz=0; iz<az.n-1; iz++) {
	for (jj=0; jj<nru[iz]; jj++) {
	    smu[iz][jj] = 0.5*(smu[iz][jj]+smu[iz+1][jj]);
	}
    }

    /*------------------------------------------------------------*/
    /* reflectivity */
    rr = sf_floatalloc2(ax.n,ay.n);

    /* allocate wavefield storage */
    ud = sf_complexalloc2(ax.n,ay.n);
    uu = sf_complexalloc2(ax.n,ay.n);
}

/*------------------------------------------------------------*/

void srmod_close(void)
/*< free allocated storage >*/
{
    ssr_close();
    
    free( *ud); free( ud);
    free( *uu); free( uu);

    free( *ss); free( ss);
    free( *so); free( so);

    free( *smd); free( smd);
    ;            free( nrd);
    free( *smu); free( smu);
    ;            free( nru);

    free( *rr); free( rr);
}

/*------------------------------------------------------------*/

void srmod_aloc()
/*< allocate migration storage >*/
{
    qq = sf_floatalloc2(ax.n,ay.n);
}

void srmod_free()
/*< free migration storage >*/
{
    free( *qq); free( qq);
}

/*------------------------------------------------------------*/
void srmod(fslice dwfl /* source   data [nw][ny][nx] */,
	   fslice uwfl /* receiver data [nw][ny][nx] */,
	   fslice refl,
	   fslice wfld
    )
/*< Apply S/R modeling >*/
{
    int iz,iw,ix,iy,ilx,ily;
    float complex w;
    
    for (iw=0; iw<aw.n; iw++) {
	if(verb) sf_warning("iw=%3d of %3d",iw+1,aw.n);
	w = eps*aw.d + I*(aw.o+iw*aw.d);

	/* downgoing wavefield */
	fslice_get(dwfl,iw,ud[0]); taper2(ud);

	fslice_put(wfld,0,ud[0]);

	fslice_get(slowd,0,so[0]);
	for (iz=0; iz<az.n-1; iz++) {
	    fslice_get(slowd,iz+1,ss[0]);

	    ssr_ssf(w,ud,so,ss,nrd[iz],smd[iz]);
	    SOOP( so[ily][ilx] = ss[ily][ilx]; );

	    fslice_put(wfld,iz+1,ud[0]);
	}
	
	/* upgoing wavefield */
	LOOP( uu[iy][ix] = 0; );

	fslice_get(slowu,az.n-1,so[0]);
	for (iz=az.n-1; iz>0; iz--) {
	    fslice_get(slowu,iz-1,ss[0]);

	    fslice_get(wfld,iz,ud[0]); 
	    fslice_get(refl,iz,rr[0]); /* reflectivity */
	    LOOP( ud[iy][ix] *= rr[iy][ix]; );
	    LOOP( uu[iy][ix] += ud[iy][ix]; );

	    ssr_ssf(w,uu,so,ss,nru[iz],smu[iz]);
	    SOOP( so[ily][ilx] = ss[ily][ilx]; );
	}
	fslice_put(uwfl,iw,uu[0]);
    } /* iw */
}
