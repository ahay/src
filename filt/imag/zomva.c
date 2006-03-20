/* 3-D SSR MVA
 * pcs 2005
 */

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

#include <rsf.h>
/*^*/

#include "zomva.h"

#include "zomig.h"
#include "taper.h"
#include "slowref.h"

#include "ssr.h"
#include "lsr.h"

#include "slice.h"
/*^*/

#define LOOP(a) for(imy=0;imy<amy.n;imy++){ for(imx=0;imx<amx.n;imx++){ {a} }}
#define SOOP(a) for(ily=0;ily<aly.n;ily++){ for(ilx=0;ilx<alx.n;ilx++){ {a} }}

static sf_axa amz,aw,alx,aly,amx,amy;
static bool verb;
static float eps;
static float twoway;

static float         **sm; /* reference slowness squared */
static float         **ss; /* slowness */
static float         **so; /* slowness */
static int            *nr; /* number of references */
static fslice       Bslow; /* slowness slice */
static fslice       Bwfld; /* wavefield slice */

static float complex **bw; /* wavefield */

static float complex **dw; /* wavefield */
static float complex **pw;
static float complex **pwsum;

static float complex **ds; /* slowness */
static float complex **ps;
static float complex **pssum;

void zomva_init(bool verb_,
		float eps_,
		bool twoway_,
		float dtmax,
		sf_axis aw_      /* frequency */,
		sf_axis amx_     /* i-line (data/image) */,
		sf_axis amy_     /* x-line (data/image) */,
		sf_axis amz_     /* depth */,
		sf_axis alx_     /* i-line (slowness) */,
		sf_axis aly_     /* x-line (slowness) */,
		int tmx, int tmy /* taper size */,
		int pmx, int pmy /* padding in the k domain */,
		int nrmax        /* maximum number of references */,
		fslice slow_,
		fslice wfld_)
/*< initialize >*/
{

    int   ilx, ily, iz, jj;
    float dsmax;

    verb=verb_;
    eps = eps_;

    aw = sf_nod(aw_);
    amx= sf_nod(amx_);
    amy= sf_nod(amy_);
    amz= sf_nod(amz_);
    alx= sf_nod(alx_);
    aly= sf_nod(aly_);

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

    /* LSR */
    lsr_init(amz_ ,
	     amx_,amy_,
	     pmx ,pmy);

    /* precompute taper */
    taper2_init(amy.n,
		amx.n ,
		SF_MIN(tmy,amy.n-1),
		SF_MIN(tmx,amx.n-1), 
		true,
		true);
    
    /* compute reference slowness */
    if(twoway_) { twoway = 2;
    } else {      twoway = 1;
    }

    Bwfld = wfld_;
    Bslow = slow_;

    ss = sf_floatalloc2(alx.n,aly.n); /* slowness */
    so = sf_floatalloc2(alx.n,aly.n); /* slowness */
    sm = sf_floatalloc2 (nrmax,amz.n); /* ref slowness squared*/
    nr = sf_intalloc          (amz.n); /* nr of ref slownesses */
    for (iz=0; iz<amz.n; iz++) {
	fslice_get(Bslow,iz,ss[0]);
	SOOP( ss[ily][ilx] *= twoway; ); /* 2-way time */

	nr[iz] = slowref(nrmax,dsmax,alx.n*aly.n,ss[0],sm[iz]);
	if (verb) sf_warning("nr[%d]=%d",iz,nr[iz]);
    }
    for (iz=0; iz<amz.n-1; iz++) {
	for (jj=0; jj<nr[iz]; jj++) {
	    sm[iz][jj] = 0.5*(sm[iz][jj]+sm[iz+1][jj]);
	}
    }

    /* allocate wavefield storage */
    bw = sf_complexalloc2(amx.n,amy.n);

    ds = sf_complexalloc2(amx.n,amy.n);
    ps = sf_complexalloc2(amx.n,amy.n);

    dw = sf_complexalloc2(amx.n,amy.n);
    pw = sf_complexalloc2(amx.n,amy.n);

    pwsum = sf_complexalloc2(amx.n,amy.n);
    pssum = sf_complexalloc2(amx.n,amy.n);
}

/*------------------------------------------------------------*/

void zomva_close(void)
/*< free allocated storage >*/
{
    ssr_close();
    lsr_close();

    free( *bw); free( bw);

    free( *ds); free( ds);
    free( *ps); free( ps);

    free( *dw); free( dw);
    free( *pw); free( pw);

    free( *pwsum); free( pwsum);
    free( *pssum); free( pssum);

    free( *ss); free( ss);
    free( *so); free( so);
    free( *sm); free( sm);
    ;           free( nr);
}

/*------------------------------------------------------------*/

void zomva(bool inv     /* forward/adjoint flag */, 
	   fslice Pslow /* slowness perturbation [nz][nmy][nmx] */,
	   fslice Pimag /*    image perturbation [nz][nmy][nmx] */)
/*< Apply forward/adjoint ZO MVA >*/
{
    int iz,iw,imy,imx,ilx,ily;
    float complex w;

    if(inv) {
	LOOP( ps[imy][imx] = 0.0; );
	for (iz=0; iz<amz.n; iz++) {
	    fslice_put(Pslow,iz,ps[0]);
	}
    } else {
	LOOP( pwsum[imy][imx] = 0.0; );
	for (iz=0; iz<amz.n; iz++) {
	    fslice_put(Pimag,iz,pwsum[0]);
	}
    }
    
    /* loop over frequencies w */
    for (iw=0; iw<aw.n; iw++) {
	if (verb) sf_warning ("iw=%3d of %3d",iw+1,aw.n);

	LOOP( dw[imy][imx]=0.; );
	
	if (inv) { /* adjoint: image -> slowness */
	    w = eps*aw.d + I*(aw.o+iw*aw.d); /* causal */

	    for (iz=amz.n-1; iz>0; iz--) {
		/* background */
		fslice_get(Bslow,iz,so[0]);	    
		SOOP( so[ily][ilx] *= twoway; ); /* 2-way time */
		fslice_get(Bwfld,iw*amz.n+iz,bw[0]);

		if(iz>0) { /* continuation */
		    fslice_get(Bslow,iz-1,ss[0]);
		    SOOP( ss[ily][ilx] *= twoway; );
		    ssr_ssf(w,dw,so,ss,nr[iz-1],sm[iz-1]);
		    SOOP( so[ily][ilx] = ss[ily][ilx]; );
		} /* end continuation */

		/* scattering dI -> dS */
		fslice_get(Pimag,iz,pwsum[0]);

		LOOP(dw[imy][imx] += pwsum[imy][imx]; );

		lsr_w2s(w,bw,so,dw,ps);

		fslice_get(Pslow,iz,pssum[0]);
		LOOP(pssum[imy][imx] += ps[imy][imx];)
		fslice_put(Pslow,iz,pssum[0]);
		/* end scattering */
	    }

	} else {   /* forward: slowness -> image */
	    w = eps*aw.d - I*(aw.o+iw*aw.d); /* anti-causal */

	    for (iz=0; iz<amz.n-1; iz++) {
		/* background */
		fslice_get(Bslow,iz,so[0]);	    
		SOOP( so[ily][ilx] *= twoway; ); /* 2-way time */
		fslice_get(Bwfld,iw*amz.n+iz,bw[0]);

		/* scattering dS -> dI */
		fslice_get(Pslow,iz,ps[0]);

		lsr_s2w(w,bw,so,pw,ps);

		LOOP(dw[imy][imx] += pw[imy][imx]; );

		fslice_get(Pimag,iz,pwsum[0]);
		LOOP(pwsum[imy][imx] += dw[imy][imx]; );
		fslice_put(Pimag,iz,pwsum[0]);
		/* end scattering */

		if(iz<amz.n-1) { /* continuation */
		    fslice_get(Bslow,iz+1,ss[0]);
		    SOOP( ss[ily][ilx] *= twoway; );
		    ssr_ssf(w,dw,so,ss,nr[iz],sm[iz]);
		    SOOP( so[ily][ilx] = ss[ily][ilx]; );
		} /* end continuation */
	    }
	}
    }
}
