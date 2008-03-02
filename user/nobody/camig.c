/* 3-D common-azimuth migration/modeling using extended split-step */
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

#include "camig.h"
#include "taper.h"
#include "slowref.h"
#include "cam.h"

#include "slice.h"
/*^*/

#define LOOP(a) for(ihx=0;ihx<ahx.n;ihx++){ \
                for(imy=0;imy<amy.n;imy++){ \
                for(imx=0;imx<amx.n;imx++){ {a} }}}
#define SOOP(a) for(ily=0;ily<aly.n;ily++){ \
                for(ilx=0;ilx<alx.n;ilx++){ {a} }}

static sf_axa az,aw,alx,aly,amx,amy,ahx,ae;
static bool verb;
static float eps;

static float         ***qq; /* image */
static sf_complex ***wx; /* wavefield x */
static float          **sm; /* reference slowness squared */
static float          **ss; /* slowness */
static float          **so; /* slowness */
static int             *nr; /* number of references */
static fslice         slow; /* slowness slice */

/*------------------------------------------------------------*/

void camig_init(bool verb_,
		float eps_,
		float dtmax,
		sf_axis az_                /* depth */,
		sf_axis aw_                /* frequency */,
		sf_axis ae_                /* experiment */,
		sf_axis amx_               /* i-line (data/image) */,
		sf_axis amy_               /* x-line (data/image) */,
		sf_axis ahx_               /* half-offset */,
		sf_axis alx_               /* i-line (slowness) */,
		sf_axis aly_               /* x-line (slowness) */,
		int tmx, int tmy, int thx /* taper size */,
		int pmx, int pmy, int phx /* padding in the k domain */,
		int nrmax                 /* maximum number of references */,
		fslice slow_)
/*< initialize >*/
{
    int   iz, j;
    float dsmax;

    verb=verb_;
    eps = eps_;

    az = sf_nod(az_);
    aw = sf_nod(aw_);
    ae = sf_nod(ae_);
    amx= sf_nod(amx_);
    amy= sf_nod(amy_);
    alx= sf_nod(alx_);
    aly= sf_nod(aly_);
    ahx= sf_nod(ahx_);
    
    /* from hertz to radian */
    aw.d *= 2.*SF_PI; 
    aw.o *= 2.*SF_PI;

    dsmax  = dtmax/az.d;

    /* CAM */
    cam_init(az_ ,
	     amx_,amy_,ahx_,
	     alx_,aly_,
	     pmx ,pmy, phx,
	     tmx ,tmy, thx,
	     dsmax);

    /* precompute taper */
    taper3_init(ahx.n,
		amy.n,
		amx.n,
		SF_MIN(thx,ahx.n-1),
		SF_MIN(tmy,amy.n-1),
		SF_MIN(tmx,amx.n-1),
		false,true,true);
		
    /* compute reference slowness */
    ss = sf_floatalloc2(alx.n,aly.n);  /* slowness */
    so = sf_floatalloc2(alx.n,aly.n);  /* slowness */
    sm = sf_floatalloc2 (nrmax,az.n);  /* reference slowness squared*/
    nr = sf_intalloc          (az.n);  /* number of reference slownesses */
    slow = slow_;
    for (iz=0; iz<az.n; iz++) {
	fslice_get(slow,iz,ss[0]);
	
	nr[iz] = slowref(nrmax,dsmax,alx.n*aly.n,ss[0],sm[iz]);
	if (verb) sf_warning("nr[%d]=%d",iz,nr[iz]);
    }
    for (iz=0; iz<az.n-1; iz++) {
	for (j=0; j<nr[iz]; j++) {
	    sm[iz][j] = 0.5*(sm[iz][j]+sm[iz+1][j]);
	}
    }

    /* allocate wavefield storage */
    wx = sf_complexalloc3(amx.n,amy.n,ahx.n);
}

/*------------------------------------------------------------*/

void camig_close(void)
/*< free allocated storage >*/
{
    cam_close();

    free(**wx); free( *wx); free( wx);
    ;           free( *ss); free( ss);
    ;           free( *so); free( so);
    ;           free( *sm); free( sm);
    ;           ;           free( nr);
}

/*------------------------------------------------------------*/

void camig_aloc()
/*< allocate migration storage >*/
{
    qq = sf_floatalloc3(amx.n,amy.n,ahx.n);
}

void camig_free()
/*< free migration storage >*/
{
    free(**qq); free( *qq); free( qq);
}
/*------------------------------------------------------------*/

void camig(bool inv  /* forward/adjoint flag */, 
	   fslice data /* data  [nw][nhx][nmy][nmx] */,
	   fslice imag /* image [nz][nhx][nmy][nmx] */)
/*< Apply migration/modeling >*/
{
    int iz,iw,imy,imx,ihx,ilx,ily;
    sf_complex w;

    if (!inv) { /* prepare image for migration */
	LOOP( qq[ihx][imy][imx] = 0.0; );
	for (iz=0; iz<az.n; iz++) {
	    fslice_put(imag,iz,qq[0][0]);
	}
    }
    
    /* loop over frequencies w */
    for (iw=0; iw<aw.n; iw++) {
/*	if (verb) sf_warning ("iw=%3d of %3d",iw+1,aw.n);*/

	if (inv) { /* MODELING */
	    w = sf_cmplx(eps*aw.d,
			 aw.o+iw*aw.d); /* +1 for upward continuation */
	    LOOP( wx[ihx][imy][imx] = sf_cmplx(0,0); );  

	    fslice_get(slow,az.n-1,so[0]);
	    for (iz=az.n-1; iz>0; iz--) {
		if (verb) sf_warning ("iw=%3d of %3d:   iz=%3d of %3d",
				      iw+1,aw.n,iz+1,az.n);
		fslice_get(imag,iz,qq[0][0]);
#ifdef SF_HAS_COMPLEX_H
		LOOP( wx[ihx][imy][imx] +=
		      qq[ihx][imy][imx]; );
#else
		LOOP( wx[ihx][imy][imx].r += 
		      qq[ihx][imy][imx]; );
#endif		
		/* upward continuation */
		fslice_get(slow,iz-1,ss[0]);
		cam_ssf(w,wx,so,ss,nr[iz-1],sm[iz-1]);
		SOOP( so[ily][ilx] = ss[ily][ilx]; );
	    }

	    fslice_get(imag,0,qq[0][0]);      /*     imaging @ iz=0 */
#ifdef SF_HAS_COMPLEX_H
	    LOOP( wx[ihx][imy][imx]  += 
		  qq[ihx][imy][imx]; );	
#else
	    LOOP( wx[ihx][imy][imx].r  += 
		  qq[ihx][imy][imx]; );
#endif	
	    
	    taper3(wx);
	    fslice_put(data,iw,wx[0][0]);    /* output data @ iz=0 */
	    
	} else { /* MIGRATION */
	    w = sf_cmplx(eps*aw.d,
			 -(aw.o+iw*aw.d)); /* -1 for downward continuation */

	    fslice_get(data,iw,wx[0][0]);    /*  input data @ iz=0 */
	    taper3(wx);

	    fslice_get(imag,0,qq[0][0]);      /*     imaging @ iz=0 */
	    LOOP(;      qq[ihx][imy][imx] += 
		 crealf(wx[ihx][imy][imx] ); );
	    fslice_put(imag,0,qq[0][0]);

	    fslice_get(slow,0,so[0]);
	    for (iz=0; iz<az.n-1; iz++) {
		if (verb) sf_warning ("iw=%3d of %3d:   iz=%3d of %3d",iw+1,aw.n,iz+1,az.n-1);

		/* downward continuation */
		fslice_get(slow,iz+1,ss[0]);
		cam_ssf(w,wx,so,ss,nr[iz],sm[iz]);
		SOOP( so[ily][ilx] = ss[ily][ilx]; );

		fslice_get(imag,iz+1,qq[0][0]); /* imaging */
		LOOP(;      qq[ihx][imy][imx] += 
		     crealf(wx[ihx][imy][imx] ); );
		fslice_put(imag,iz+1,qq[0][0]);
	    }

	} /* else */
    } /* iw */
}

/*------------------------------------------------------------*/

void cadtm(bool inv    /* forward/adjoint flag */, 
	   fslice data /* data [nw][nhx][nmy][nmx] */,
	   fslice wfld /* wfld [nw][nhx][nmy][nmx] */)
/*< Apply upward/downward datuming >*/
{
    int iz,iw,ie, ilx,ily;
    sf_complex w;
    
    for(ie=0; ie<ae.n; ie++) {
	
	/* loop over frequencies w */
	for (iw=0; iw<aw.n; iw++) {
	    if (verb) sf_warning ("iw=%3d of %3d:   ie=%3d of %3d",
				  iw+1,aw.n,ie+1,ae.n);
	    
	    if (inv) { /* UPWARD DATUMING */
		w = sf_cmplx(eps*aw.d,aw.o+iw*aw.d);
		
		fslice_get(wfld,iw+ie*aw.n,wx[0][0]);
		taper3(wx);
		
		fslice_get(slow,az.n-1,so[0]);
		for (iz=az.n-1; iz>0; iz--) {
		    fslice_get(slow,iz-1,ss[0]);
		    cam_ssf(w,wx,so,ss,nr[iz-1],sm[iz-1]);
		    SOOP( so[ily][ilx] = ss[ily][ilx]; );
		}
		
		taper3(wx);
		fslice_put(data,iw+ie*aw.n,wx[0][0]);
	    } else { /* DOWNWARD DATUMING */
		w = sf_cmplx(eps*aw.d,-(aw.o+iw*aw.d));
		
		fslice_get(data,iw+ie*aw.n,wx[0][0]);
		taper3(wx);
		
		fslice_get(slow,0,so[0]);
		for (iz=0; iz<az.n-1; iz++) {
		    fslice_get(slow,iz+1,ss[0]);
		    cam_ssf(w,wx,so,ss,nr[iz],sm[iz]);
		    SOOP( so[ily][ilx] = ss[ily][ilx]; );
		}
		
		taper3(wx);
		fslice_put(wfld,iw+ie*aw.n,wx[0][0]);
	    } /* else */
	} /* iw */
    } /* ie */
}

/*------------------------------------------------------------*/

void cawfl(fslice data /*      data [nw][nhx][nmy][nmx] */,
	   fslice wfld /* wavefield [nw][nhx][nmy][nmx] */)
/*< Save wavefield from downward continuation >*/
{
    int iz,iw, ilx,ily;
    sf_complex w;

    /* loop over frequencies w */
    for (iw=0; iw<aw.n; iw++) {
/*	if (verb) sf_warning ("iw=%3d of %3d",iw+1,aw.n);*/
	w = sf_cmplx(eps*aw.d,aw.o+iw*aw.d);

	fslice_get(data,iw,wx[0][0]);
	taper3(wx);
	
	taper3(wx);
	fslice_put(wfld,iw*az.n,wx[0][0]);

	fslice_get(slow,0,so[0]);
	for (iz=0; iz<az.n-1; iz++) {	
	    if (verb) sf_warning ("iw=%3d of %3d:   iz=%3d of %3d",
				  iw+1,aw.n,iz+1,az.n);

	    fslice_get(slow,iz+1,ss[0]);
	    cam_ssf(w,wx,so,ss,nr[iz],sm[iz]);
	    SOOP( so[ily][ilx] = ss[ily][ilx]; );

	    taper3(wx);
	    fslice_put(wfld,iw*az.n+iz+1,wx[0][0]);
	}	    
    } /* iw */
}

/*------------------------------------------------------------*/

