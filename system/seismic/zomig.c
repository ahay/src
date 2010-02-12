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

#include "zomig.h"
#include "taper.h"
#include "slowref.h"
#include "ssr.h"

#define LOOP(a) for(imy=0;imy<amy.n;imy++){ for(imx=0;imx<amx.n;imx++){ {a} }}
#define SOOP(a) for(ily=0;ily<aly.n;ily++){ for(ilx=0;ilx<alx.n;ilx++){ {a} }}

static sf_axa az,aw,alx,aly,amx,amy,ae;
static bool verb, incore;
static float eps;
static float twoway;

static float         **qq; /* image */
static sf_complex **wx; /* wavefield x */
static float         **sm; /* reference slowness squared */
static float         **ss; /* slowness */
static float         **so; /* slowness */
static int            *nr; /* number of references */
static sf_fslice        slow; /* slowness slice */

static float       ***rrr;
static float       ***sss;

/*------------------------------------------------------------*/

void zomig_init(bool verb_,
		bool incore_,  /* keep shot wavefield in core */
		float eps_,
		bool twoway_,
		float dtmax,
		sf_axis az_      /* depth */,
		sf_axis aw_      /* frequency */,
		sf_axis ae_      /* experiment */,
		sf_axis amx_     /* i-line (data/image) */,
		sf_axis amy_     /* x-line (data/image) */,
		sf_axis alx_     /* i-line (slowness) */,
		sf_axis aly_     /* x-line (slowness) */,
		int tmx, int tmy /* taper size */,
		int pmx, int pmy /* padding in the k domain */,
		int nrmax        /* maximum number of references */,
		sf_fslice slow_)
/*< initialize >*/
{
    int   ilx, ily, iz, jj;
    float dsmax;

    verb   = verb_;
    incore = incore_;
    eps    = eps_;

    az = sf_nod(az_);
    aw = sf_nod(aw_);
    ae = sf_nod(ae_);
    amx= sf_nod(amx_);
    amy= sf_nod(amy_);
    alx= sf_nod(alx_);
    aly= sf_nod(aly_);

    /* from hertz to radian */
    aw.d *= 2.*SF_PI;
    aw.o *= 2.*SF_PI;

    dsmax  = dtmax/az.d;

    /* SSR */
    ssr_init(az_ ,
	     amx_,amy_,
	     alx_,aly_,
	     pmx ,pmy,
	     tmx ,tmy,
	     dsmax);

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

    ss = sf_floatalloc2(alx.n,aly.n); /* slowness */
    so = sf_floatalloc2(alx.n,aly.n); /* slowness */
    sm = sf_floatalloc2 (nrmax,az.n); /* ref slowness squared*/
    nr = sf_intalloc          (az.n); /* nr of ref slownesses */
    slow = slow_;
    for (iz=0; iz<az.n; iz++) {
	sf_fslice_get(slow,iz,ss[0]);
	SOOP( ss[ily][ilx] *= twoway; ); /* 2-way time */

	nr[iz] = slowref(nrmax,dsmax,alx.n*aly.n,ss[0],sm[iz]);
	if(verb) sf_warning("nr[%d]=%d",iz,nr[iz]);
    }
    for (iz=0; iz<az.n-1; iz++) {
	for (jj=0; jj<nr[iz]; jj++) {
	    sm[iz][jj] = 0.5*(sm[iz][jj]+sm[iz+1][jj]);
	    if(verb) sf_warning("s[%d][%d]=%g",iz,jj,sm[iz][jj]);
	}
    }
    
    if(incore) {
	sss = sf_floatalloc3  (alx.n,aly.n,az.n);
    }

    /* allocate wavefield storage */
    wx = sf_complexalloc2(amx.n,amy.n);
}

/*------------------------------------------------------------*/

void zomig_close(void)
/*< free allocated storage >*/
{
    ssr_close();

    free( *wx); free( wx);
    free( *ss); free( ss);
    free( *so); free( so);
    free( *sm); free( sm);
    ;           free( nr);

    if(incore) {
	free(**sss); free(*sss); free(sss);
    }
}

/*------------------------------------------------------------*/

void zomig_aloc()
/*< allocate migration storage >*/
{
    qq = sf_floatalloc2(amx.n,amy.n);
    if(incore) {
	rrr = sf_floatalloc3(amx.n,amy.n,az.n);
    }
}

void zomig_free()
/*< free migration storage >*/
{
    free( *qq); free( qq);
    if(incore) {
	free(**rrr); free( *rrr); free( rrr);
    }
}

/*------------------------------------------------------------*/

void zomig(bool inv    /* forward/adjoint flag */, 
	   sf_fslice data /* data  [nw][nmy][nmx] */,
	   sf_fslice imag /* image [nz][nmy][nmx] */)
/*< Apply migration/modeling >*/
{
    int iz,iw,imy,imx,ilx,ily;
    sf_complex w;

    if(!inv) { /* prepare image for migration */
	LOOP( qq[imy][imx] = 0.0; );
	for (iz=0; iz<az.n; iz++) {
	    sf_fslice_put(imag,iz,qq[0]);
	}
    }

    if(incore) {
	for( iz=0; iz<az.n; iz++) {
	    LOOP( rrr[iz][imy][imx] = 0;);

	    sf_fslice_get(slow,iz,so[0]);
	    SOOP( sss[iz][ily][ilx] = so[ily][ilx]*twoway;);
	}

	if(inv) { /* MODELING */
	    
	    for (iz=0; iz<az.n; iz++) {
		sf_fslice_get(imag,iz,qq[0]);
		LOOP( rrr[iz][imy][imx] = qq[imy][imx]; );
	    }
	    for (iw=0; iw<aw.n; iw++) { /* frequency loop */
		if(verb) sf_warning ("iw=%3d of %3d (in core)",iw+1,aw.n);

		w = sf_cmplx(eps*aw.d,aw.o+iw*aw.d); /* causal */
		
		LOOP( wx[imy][imx] = sf_cmplx(0.,0.); );		
		for (iz=az.n-1; iz>0; iz--) {
                    /* imaging */
#ifdef SF_HAS_COMPLEX_H
		    LOOP( wx[imy][imx] += rrr[iz][imy][imx]; );    
#else
		    LOOP( wx[imy][imx].r += rrr[iz][imy][imx]; );
#endif    
		    /* extrapolation */
		    ssr_ssf(w,wx,sss[iz],sss[iz-1],nr[iz-1],sm[iz-1]); 
		}
#ifdef SF_HAS_COMPLEX_H		
		LOOP( wx[imy][imx] += rrr[0][imy][imx]; );
#else
		LOOP( wx[imy][imx].r += rrr[0][imy][imx]; );
#endif
		taper2(wx);
		sf_fslice_put(data,iw,wx[0]);
	    }
	} else { /* MIGRATION */

	    for (iw=0; iw<aw.n; iw++) { /* frequency loop */
		if(verb) sf_warning ("iw=%3d of %3d (in core)",iw+1,aw.n);
		w = sf_cmplx(eps*aw.d,-(aw.o+iw*aw.d)); /* anti-causal */
		
		sf_fslice_get(data,iw,wx[0]); 
		taper2(wx);
		LOOP( rrr[0][imy][imx] += crealf(wx[imy][imx]); );
		
		for (iz=0; iz<az.n-1; iz++) {		    
                   /* extrapolation */
		    ssr_ssf(w,wx,sss[iz],sss[iz+1],nr[iz],sm[iz]);        
		    /* imaging */
		    LOOP( rrr[iz+1][imy][imx] += crealf(wx[imy][imx]); ); 
		}
	    }
	    for (iz=0; iz<az.n; iz++) {
		LOOP(qq[imy][imx] = rrr[iz][imy][imx]; );
		sf_fslice_put(imag,iz,qq[0]);
	    }
	} /* else inv */

    } else {
	/* loop over frequencies w */
	for (iw=0; iw<aw.n; iw++) {
	    if(verb) sf_warning ("iw=%3d of %3d",iw+1,aw.n);
	    
	    if(inv) { /* MODELING */
		w = sf_cmplx(eps*aw.d,aw.o+iw*aw.d); /* causal */
		LOOP( wx[imy][imx] = sf_cmplx(0.,0.); );  
		
		sf_fslice_get(slow,az.n-1,so[0]);
		SOOP( so[ily][ilx] *= twoway; ); /* 2-way time */
		for (iz=az.n-1; iz>0; iz--) {
		    sf_fslice_get(imag,iz,qq[0]);
#ifdef SF_HAS_COMPLEX_H
		    LOOP( wx[imy][imx] +=
			  qq[imy][imx]; );
#else
		     LOOP( wx[imy][imx].r +=
			   qq[imy][imx]; );
#endif
		    
		    /* upward continuation */
		    sf_fslice_get(slow,iz-1,ss[0]);
		    SOOP( ss[ily][ilx] *= twoway; ); /* 2-way time */
		    ssr_ssf(w,wx,so,ss,nr[iz-1],sm[iz-1]);
		    SOOP( so[ily][ilx] = ss[ily][ilx]; );
		}
		
		sf_fslice_get(imag,0,qq[0]);      /*     imaging @ iz=0 */
#ifdef SF_HAS_COMPLEX_H
		LOOP( wx[imy][imx]  += 
		      qq[imy][imx]; );	
#else
		LOOP( wx[imy][imx].r  += 
		      qq[imy][imx]; );
#endif
		
		taper2(wx);
		sf_fslice_put(data,iw,wx[0]);    /* output data @ iz=0 */
		
	    } else { /* MIGRATION */
		w = sf_cmplx(eps*aw.d,-(aw.o+iw*aw.d)); /* anti-causal */

		sf_fslice_get(data,iw,wx[0]);    /*  input data @ iz=0 */
		taper2(wx);
		
		sf_fslice_get(imag,0,qq[0]);      /*     imaging @ iz=0 */
		LOOP(;      qq[imy][imx] += 
		     crealf(wx[imy][imx] ); );
		sf_fslice_put(imag,0,qq[0]);
		
		sf_fslice_get(slow,0,so[0]);	    
		SOOP( so[ily][ilx] *= twoway; ); /* 2-way time */
		for (iz=0; iz<az.n-1; iz++) {
		    
		    /* downward continuation */
		    sf_fslice_get(slow,iz+1,ss[0]);
		    SOOP( ss[ily][ilx] *= twoway; ); /* 2-way time */
		    ssr_ssf(w,wx,so,ss,nr[iz],sm[iz]);
		    SOOP( so[ily][ilx] = ss[ily][ilx]; );
		    
		    sf_fslice_get(imag,iz+1,qq[0]); /* imaging */
		    LOOP(;      qq[imy][imx] += 
			 crealf(wx[imy][imx] ); );
		    sf_fslice_put(imag,iz+1,qq[0]);
		} /* iz */
	    } /* else */
	} /* iw */
    } /* incore */
}

/*------------------------------------------------------------*/

void zodtm(bool inv    /* forward/adjoint flag */, 
	   bool causal,
	   sf_fslice data /* data [nw][nmy][nmx] */,
	   sf_fslice wfld /* wfld [nw][nmy][nmx] */)
/*< Apply upward/downward datuming >*/
{
    int iz,iw,ie, ilx,ily;
    sf_complex w;

    for(ie=0; ie<ae.n; ie++) {
	
	/* loop over frequencies w */
	for (iw=0; iw<aw.n; iw++) {
	    if(verb) sf_warning ("iw=%3d of %3d:   ie=%3d of %3d",
				 iw+1,aw.n,ie+1,ae.n);
	    
	    if(inv) { /* UPWARD DATUMING */
		if(causal) {
		    w = sf_cmplx(eps*aw.d,-(aw.o+iw*aw.d));
		} else {
		    w = sf_cmplx(eps*aw.d,  aw.o+iw*aw.d );
		}

		sf_fslice_get(wfld,iw+ie*aw.n,wx[0]);
		taper2(wx);
		
		sf_fslice_get(slow,az.n-1,so[0]);
		SOOP( so[ily][ilx] *= twoway; );     /* 2-way time */
		for (iz=az.n-1; iz>0; iz--) {
		    sf_fslice_get(slow,iz-1,ss[0]);
		    SOOP( ss[ily][ilx] *= twoway; ); /* 2-way time */
		    ssr_ssf(w,wx,so,ss,nr[iz-1],sm[iz-1]);
		    SOOP( so[ily][ilx] = ss[ily][ilx]; );
		}
		
		taper2(wx);
		sf_fslice_put(data,iw+ie*aw.n,wx[0]);
	    } else { /* DOWNWARD DATUMING */
		if(causal) {
		    w = sf_cmplx(eps*aw.d,  aw.o+iw*aw.d );
		} else {
		    w = sf_cmplx(eps*aw.d,-(aw.o+iw*aw.d));
		}

		sf_fslice_get(data,iw+ie*aw.n,wx[0]);
		taper2(wx);
		
		sf_fslice_get(slow,0,so[0]);
		SOOP( so[ily][ilx] *= twoway; );     /* 2-way time */
		for (iz=0; iz<az.n-1; iz++) {
		    sf_fslice_get(slow,iz+1,ss[0]);
		    SOOP( ss[ily][ilx] *= twoway; ); /* 2-way time */
		    ssr_ssf(w,wx,so,ss,nr[iz],sm[iz]);
		    SOOP( so[ily][ilx] = ss[ily][ilx]; );
		}
		
		taper2(wx);
		sf_fslice_put(wfld,iw+ie*aw.n,wx[0]);
	    } /* else */
	} /* iw */
    } /* ie */
}

/*------------------------------------------------------------*/

void zowfl(bool inv    /* forward/adjoint flag */, 
	   bool causal,
	   sf_fslice data /*      data [nw][nmy][nmx] */,
	   sf_fslice wfld /* wavefield [nw][nmy][nmx] */)
/*< Save wavefield from downward continuation >*/
{
    int iz,iw, ilx,ily;
    sf_complex w;

    /* loop over frequencies w */
    for (iw=0; iw<aw.n; iw++) {
	if(verb) sf_warning ("iw=%3d of %3d",iw+1,aw.n);

	if(inv) { /*   UPWARD EXTRAPOLATION */
	    if(causal) {
		w = sf_cmplx(eps*aw.d,-(aw.o+iw*aw.d));
	    } else {
		w = sf_cmplx(eps*aw.d,  aw.o+iw*aw.d );
	    }
	    
	    sf_fslice_get(data,iw,wx[0]);
	    taper2(wx);

	    taper2(wx);
	    sf_fslice_put(wfld,iw*az.n+az.n-1,wx[0]);

	    sf_fslice_get(slow,az.n-1,so[0]);
	    SOOP( so[ily][ilx] *= twoway; );     /* 2-way time */
	    for (iz=az.n-1; iz>0; iz--) {
		sf_fslice_get(slow,iz-1,ss[0]);
		SOOP( ss[ily][ilx] *= twoway; ); /* 2-way time */
		ssr_ssf(w,wx,so,ss,nr[iz-1],sm[iz-1]);
		SOOP( so[ily][ilx] = ss[ily][ilx]; );

		taper2(wx);
		sf_fslice_put(wfld,iw*az.n+iz-1,wx[0]);
	    }

	} else {  /* DOWNWARD EXTRAPOLATION */
	    if(causal) {
		w = sf_cmplx(eps*aw.d,  aw.o+iw*aw.d );
	    } else {
		w = sf_cmplx(eps*aw.d,-(aw.o+iw*aw.d));
	    }

	    sf_fslice_get(data,iw,wx[0]);
	    taper2(wx);

	    sf_fslice_put(wfld,iw*az.n,wx[0]);

	    sf_fslice_get(slow,0,so[0]);
	    SOOP( so[ily][ilx] *= twoway; );     /* 2-way time */
	    for (iz=0; iz<az.n-1; iz++) {	    
		sf_fslice_get(slow,iz+1,ss[0]);
		SOOP( ss[ily][ilx] *= twoway; ); /* 2-way time */
		ssr_ssf(w,wx,so,ss,nr[iz],sm[iz]);
		SOOP( so[ily][ilx] = ss[ily][ilx]; );
		
		sf_fslice_put(wfld,iw*az.n+iz+1,wx[0]);
	    } /* iz */
	} /* else */
    } /* iw */
}
/*------------------------------------------------------------*/

