/* 3-D SR MVA using extended split-step */
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

#include "srmva.h"
#include "taper.h"
#include "slowref.h"

#include "ssr.h"
#include "lsr.h"

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
static sf_fslice       Bslow; /* slowness slice */
static sf_fslice       Bwfls; /* wavefield slice */
static sf_fslice       Bwflr; /* wavefield slice */

static sf_complex **bw_s,**bw_r; /* wavefield */

static sf_complex **pw_s,**pw_r; /* */ 
static sf_complex **dw_s,**dw_r;
static sf_complex **pwsum;

static sf_complex **ps;
static sf_complex **ps_s,**ps_r;
static sf_complex **ds_s,**ds_r; /* slowness */
static sf_complex **pssum;

void srmva_init(bool verb_,
		float eps_,
		bool twoway_,
		float dtmax,
		sf_axis aw_      /* frequency */,
		sf_axis amx_     /* i-line (data) */,
		sf_axis amy_     /* x-line (data) */,
		sf_axis amz_     /* depth */,
		sf_axis alx_     /* i-line (slowness/image) */,
		sf_axis aly_     /* x-line (slowness/image) */,
		int tmx, int tmy /* taper size */,
		int pmx, int pmy /* padding in the k domain */,
		int    nrmax,    /* maximum number of references */
		sf_fslice slow_,
		sf_fslice wfls_,
		sf_fslice wflr_
    )
/*< initialize >*/
{

    int   iz, jj;
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

    Bwfls = wfls_;
    Bwflr = wflr_;
    Bslow = slow_;

    ss = sf_floatalloc2(alx.n,aly.n); /* slowness */
    so = sf_floatalloc2(alx.n,aly.n); /* slowness */
    sm = sf_floatalloc2 (nrmax,amz.n); /* ref slowness squared*/
    nr = sf_intalloc          (amz.n); /* nr of ref slownesses */

    /* find reference slownesses */
    for (iz=0; iz<amz.n; iz++) {
	sf_fslice_get(Bslow,iz,ss[0]);

	nr[iz] = slowref(nrmax,dsmax,alx.n*aly.n,ss[0],sm[iz]);
	if (verb) sf_warning("nr[%d]=%d",iz,nr[iz]);
    }
    for (iz=0; iz<amz.n-1; iz++) {
	for (jj=0; jj<nr[iz]; jj++) {
	    sm[iz][jj] = 0.5*(sm[iz][jj]+sm[iz+1][jj]);
	}
    }

    /* allocate wavefield storage */

    /* background */
    bw_s  = sf_complexalloc2(amx.n,amy.n);
    bw_r  = sf_complexalloc2(amx.n,amy.n);

    /* slowness perturbations */
    ds_s  = sf_complexalloc2(amx.n,amy.n);
    ds_r  = sf_complexalloc2(amx.n,amy.n);
    ps_s  = sf_complexalloc2(amx.n,amy.n);
    ps_r  = sf_complexalloc2(amx.n,amy.n);

    /* wavefield perturbations */
    pw_s  = sf_complexalloc2(amx.n,amy.n);
    pw_r  = sf_complexalloc2(amx.n,amy.n);
    dw_s  = sf_complexalloc2(amx.n,amy.n);
    dw_r  = sf_complexalloc2(amx.n,amy.n);

    ps    = sf_complexalloc2(amx.n,amy.n);
    pwsum = sf_complexalloc2(amx.n,amy.n);
    pssum = sf_complexalloc2(amx.n,amy.n);
}
/*------------------------------------------------------------*/

void srmva_close(void)
/*< free allocated storage >*/
{
    ssr_close();
    lsr_close();
    
    free( *bw_s); free( bw_s);
    free( *bw_r); free( bw_r);

    free( *ps  ); free( ps  );

    free( *ds_s); free( ds_s);
    free( *ds_r); free( ds_r);
    free( *ps_s); free( ps_s);
    free( *ps_r); free( ps_r);

    free( *pw_s); free( pw_s);
    free( *pw_r); free( pw_r);
    free( *dw_s); free( dw_s);
    free( *dw_r); free( dw_r);

    free( *pwsum); free( pwsum);
    free( *pssum); free( pssum);

    free( *ss); free( ss);
    free( *so); free( so);
    free( *sm); free( sm);
    ;           free( nr);
}

/*------------------------------------------------------------*/

void srmva(bool adj     /* forward/adjoint flag */, 
	   sf_fslice Pslow /* slowness perturbation [nz][nmy][nmx] */,
	   sf_fslice Pimag /*    image perturbation [nz][nmy][nmx] */)
/*< Apply forward/adjoint SR MVA >*/
{
    int imz,iw,imy,imx,ilx,ily;
    sf_complex ws,wr;

    if(adj) {
	LOOP( pssum[imy][imx] = sf_cmplx(0.0,0.0); );
	for (imz=0; imz<amz.n; imz++) {
	    sf_fslice_put(Pslow,imz,pssum[0]);
	}
    } else {
	LOOP( pwsum[imy][imx] = sf_cmplx(0.0,0.0); );
	for (imz=0; imz<amz.n; imz++) {
	    sf_fslice_put(Pimag,imz,pwsum[0]);
	}
    }
    
    /* loop over frequencies w */
    for (iw=0; iw<aw.n; iw++) {
	if (verb) sf_warning ("iw=%3d of %3d",iw+1,aw.n);

	LOOP( dw_s[imy][imx]=sf_cmplx(0.,0.); 
	      dw_r[imy][imx]=sf_cmplx(0.,0.); );
	
	if (adj) { /* adjoint: image -> slowness */
	    ws = sf_cmplx(eps*aw.d,-(aw.o+iw*aw.d)); /* anti-causal */
	    wr = sf_cmplx(eps*aw.d,+(aw.o+iw*aw.d)); /*      causal */

	    imz = amz.n-1;
	    sf_fslice_get(Bslow,imz,so[0]);

	    for (imz=amz.n-1; imz>0; imz--) {
		/* background */
		sf_fslice_get(Bslow,imz,so[0]);
		sf_fslice_get(Bwfls,iw*amz.n+imz,bw_s[0]);
		sf_fslice_get(Bwflr,iw*amz.n+imz,bw_r[0]);

		if(imz<amz.n-1) { /* continuation */
		    sf_fslice_get(Bslow,imz-1,ss[0]);
		    ssr_ssf(ws,dw_s,so,ss,nr[imz],sm[imz]);
		    ssr_ssf(wr,dw_r,so,ss,nr[imz],sm[imz]);
		    SOOP( so[ily][ilx] = ss[ily][ilx]; );
		}

		/* scattering dI -> dS */
		sf_fslice_get(Pimag,imz,pwsum[0]);
#ifdef SF_HAS_COMPLEX_H
		LOOP( dw_s [imy][imx] +=
		      bw_r [imy][imx] * conjf(pwsum[imy][imx]);
		      dw_r [imy][imx] += 
		      bw_s [imy][imx] *       pwsum[imy][imx];
		    );

		lsr_w2s(ws,bw_s,so,dw_s,ps_s);
		lsr_w2s(wr,bw_r,so,dw_r,ps_r);

		sf_fslice_get(Pslow,imz,pssum[0]);

		LOOP(pssum[imy][imx] += 
		     ps_s [imy][imx]+ 
		     ps_r [imy][imx];
		    );
#else
		LOOP( dw_s [imy][imx] = 
		      sf_cadd(dw_s [imy][imx],
			      sf_cmul(bw_r [imy][imx],conjf(pwsum[imy][imx])));
		      dw_r [imy][imx] =
		      sf_cadd(dw_r [imy][imx],
			      sf_cmul(bw_s [imy][imx],      pwsum[imy][imx]));
		    );

		lsr_w2s(ws,bw_s,so,dw_s,ps_s);
		lsr_w2s(wr,bw_r,so,dw_r,ps_r);

		sf_fslice_get(Pslow,imz,pssum[0]);

		LOOP(pssum[imy][imx] = 
		     sf_cadd(pssum[imy][imx], 
			     sf_cadd(ps_s [imy][imx], 
				     ps_r [imy][imx]));
		    );
#endif
		sf_fslice_put(Pslow,imz,pssum[0]);
		/* end scattering */
	    }
	} else {   /* forward: slowness -> image */
	    ws = sf_cmplx(eps*aw.d,+(aw.o+iw*aw.d)); /*      causal */
	    wr = sf_cmplx(eps*aw.d,-(aw.o+iw*aw.d)); /* anti-causal */

	    for (imz=0; imz<amz.n-1; imz++) {
		/* background */
		sf_fslice_get(Bslow,imz,so[0]);
		sf_fslice_get(Bwfls,iw*amz.n+imz,bw_s[0]);
		sf_fslice_get(Bwflr,iw*amz.n+imz,bw_r[0]);

		/* scattering dS -> dI */
		sf_fslice_get(Pslow,imz,ps[0]);

		lsr_s2w(ws,bw_s,so,pw_s,ps);
		lsr_s2w(wr,bw_r,so,pw_r,ps);

		/* cross-correlation of background and perturbation wavefields */
#ifdef SF_HAS_COMPLEX_H
		LOOP(dw_s[imy][imx] += pw_s[imy][imx]; );
		LOOP(dw_r[imy][imx] += pw_r[imy][imx]; );
		sf_fslice_get(Pimag,imz,pwsum[0]);
		LOOP(     pwsum[imy][imx] += 
		     conjf(bw_s[imy][imx]) * dw_r[imy][imx] +
		     conjf(dw_s[imy][imx]) * bw_r[imy][imx]; 
		    );
		sf_fslice_put(Pimag,imz,pwsum[0]);
#else
		LOOP(dw_s[imy][imx] = sf_cadd(dw_s[imy][imx],pw_s[imy][imx]); );
		LOOP(dw_r[imy][imx] = sf_cadd(dw_r[imy][imx],pw_r[imy][imx]); );
		sf_fslice_get(Pimag,imz,pwsum[0]);
		LOOP(     pwsum[imy][imx] =
			  sf_cadd( pwsum[imy][imx],
				   sf_cadd( sf_cmul(conjf(bw_s[imy][imx]),dw_r[imy][imx]),
					    sf_cmul(conjf(dw_s[imy][imx]),bw_r[imy][imx]))); 
		    );
		sf_fslice_put(Pimag,imz,pwsum[0]);
#endif
		/* end scattering dS -> dI */

		if(imz<amz.n-1) { /* continuation */
		    sf_fslice_get(Bslow,imz+1,ss[0]);
		    ssr_ssf(ws,dw_s,so,ss,nr[imz],sm[imz]);
		    ssr_ssf(wr,dw_r,so,ss,nr[imz],sm[imz]);
		    SOOP( so[ily][ilx] = ss[ily][ilx]; );
		} /* end if */

	    } /* loop z */
	} /* end if adjoint */
    } /* loop w*/

}
