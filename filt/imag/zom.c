/* 3-D SSR migration/modeling using extended split-step */
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

#include "zom.h"
#include "fft2.h"
#include "taper.h"
#include "slowref.h"

#include "slice.h"
/*^*/

#define LOOP(a) for(imy=0;imy<amy.n;imy++){ for(imx=0;imx<amx.n;imx++){ {a} }}
#define KOOP(a) for(imy=0;imy<bmy.n;imy++){ for(imx=0;imx<bmx.n;imx++){ {a} }}
#define SOOP(a) for(ily=0;ily<aly.n;ily++){ for(ilx=0;ilx<alx.n;ilx++){ {a} }}

#define INDEX(x,a) 0.5+(x-a.o)/a.d;
#define BOUND(i,n) (i<0) ? 0 : ( (i>n-1) ? n-1 : i );
#define  KMAP(i,n) (i<n/2.) ? (i+n/2.) : (i-n/2.);
#define X2K(a,b,p) b.n=a.n+p; b.d=2.0*SF_PI/(b.n*a.d); b.o=(1==b.n)?0:-SF_PI/a.d;

static axa az,aw,alx,aly,amx,amy;
static axa               bmx,bmy;
static bool verb;
static float eps;
static float ds, ds2;

static float         **qq; /* image */
static float complex **pk; /* wavefield k */
static float complex **wk; /* wavefield k */
static float complex **wx; /* wavefield x */
static float         **wt; /* interpolation weight */
static float         **kk; /* wavenumber  */

static float          **sm; /* reference slowness squared */
static float          **ss; /* slowness */
static float          **so; /* slowness */

static int             *jx; /* i-line index */
static int             *jy; /* x-line index */
static int             *nr; /* number of references */

static slice          slow; /* slowness slice */

/*------------------------------------------------------------*/

void zom_init(bool verb_,
	      float eps_,
	      float  dt,
	      axa az_          /* depth */,
	      axa aw_          /* frequency */,
	      axa amx_         /* i-line (data) */,
	      axa amy_         /* x-line (data) */,
	      axa alx_         /* i-line (slowness/image) */,
	      axa aly_         /* x-line (slowness/image) */,
	      int tmx, int tmy /* taper size */,
	      int pmx, int pmy /* padding in the k domain */,
	      int nrmax        /* maximum number of references */,
	      slice slow_)
/*< initialize >*/
{
    int   imy, imx, ilx, ily, iz, j;
    int   jmy, jmx;
    float  my,  mx;
    float kmy,kmx;

    slow = slow_;

    verb=verb_;
    eps = eps_;

    az = az_;
    aw = aw_;
    amx= amx_;
    amy= amy_;
    alx= alx_;
    aly= aly_;
    
    /* from hertz to radian */
    aw.d *= 2.*SF_PI; 
    aw.o *= 2.*SF_PI;

    ds  = dt/az.d;
    ds2 = ds*ds;
    ds2*= ds2;

    /* construct K-domain axes */
    X2K(amx,bmx,pmx);
    X2K(amy,bmy,pmy);
    fft2_init(bmx.n,bmy.n);

    /* allocate storage */
    ss = sf_floatalloc2   (alx.n,aly.n);  /* slowness */
    so = sf_floatalloc2   (alx.n,aly.n);  /* slowness */
    sm = sf_floatalloc2    (nrmax,az.n);  /* reference slowness squared*/
    nr = sf_intalloc             (az.n);  /* number of reference slownesses */

    wx = sf_complexalloc2 (amx.n,amy.n);  /* x wavefield */
    wt = sf_floatalloc2   (amx.n,amy.n);  /* interpolation weight */

    wk = sf_complexalloc2 (bmx.n,bmy.n);  /*        k wavefield */
    pk = sf_complexalloc2 (bmx.n,bmy.n);  /* padded k wavefield */ 

    kk= sf_floatalloc2    (bmx.n,bmy.n);  /* wavenumber */

    jy = sf_intalloc      (amy.n);  /* midpoint index */
    jx = sf_intalloc      (amx.n);  /* midpoint index */

    /* precompute indices */
    for (imy=0; imy<amy.n; imy++) {
	my = amy.o + imy*amy.d;
	ily     = INDEX( my,aly);
	jy[imy] = BOUND(ily,aly.n); /* x-line index */
    }
    for (imx=0; imx<amx.n; imx++) {
	mx = amx.o + imx*amx.d;
	ilx     = INDEX( mx,alx);
	jx[imx] = BOUND(ilx,alx.n); /* i-line index */
    }

    /* precompute wavenumbers */
    for (imy=0; imy<bmy.n; imy++) {
	jmy = KMAP(imy,bmy.n);
	kmy = bmy.o + jmy*bmy.d;
	
	for (imx=0; imx<bmx.n; imx++) {
	    jmx = KMAP(imx,bmx.n);
	    kmx = bmx.o + jmx*bmx.d;
	    
	    kk[imy][imx] = kmx*kmx+kmy*kmy;
	}
    }    

    /* precompute taper */
    taper2_init(amy.n,amx.n,
		SF_MIN(tmy,amy.n-1),
		SF_MIN(tmx,amx.n-1) );
		
    /* compute reference slowness */
    for (iz=0; iz<az.n; iz++) {
	slice_get(slow,iz,ss[0]);
	nr[iz] = slowref(nrmax,ds,alx.n*aly.n,ss[0],sm[iz]);
	if (verb) sf_warning("nr[%d]=%d",iz,nr[iz]);
    }
    for (iz=0; iz<az.n-1; iz++) {
	for (j=0; j<nr[iz]; j++) {
	    sm[iz][j] = 0.5*(sm[iz][j]+sm[iz+1][j]);
	}
    }
}

/*------------------------------------------------------------*/

void zomig_init()
/*< allocate migration storage >*/
{
    qq = sf_floatalloc2   (amx.n,amy.n);  /* image */
}

/*------------------------------------------------------------*/

void zom_close(void)
/*< free allocated storage >*/
{
    free( *pk); free( pk);
    free( *wk); free( wk);
    free( *wx); free( wx);
    free( *wt); free( wt);
    free( *ss); free( ss);
    free( *so); free( so);
    free( *sm); free( sm);
    ;           free( nr);
    free( *kk); free( kk);
    
    fft2_close();
    taper2_close();
}

/*------------------------------------------------------------*/

void zomig_close()
/*< free migration storage >*/
{
    free( *qq); free( qq);
}
/*------------------------------------------------------------*/

void zomig(bool inv  /* forward/adjoint flag */, 
	  slice data /* data  [nw][nmy][nmx] */,
	  slice imag /* image [nz][nmy][nmx] */)
/*< Apply migration/modeling >*/
{
    int iz,iw,imy,imx,ilx,ily;
    float complex w;

    if (!inv) { /* prepare image for migration */
	LOOP( qq[imy][imx] = 0.0; );
	for (iz=0; iz<az.n; iz++) {
	    slice_put(imag,iz,qq[0]);
	}
    }
    
    /* loop over frequencies w */
    for (iw=0; iw<aw.n; iw++) {
	if (verb) sf_warning ("frequency %d of %d",iw+1,aw.n);

	if (inv) { /* MODELING */
	    w = eps*aw.d + I*(aw.o+iw*aw.d); /* +1 for upward continuation */
	    LOOP( wx[imy][imx] = 0; );  

	    slice_get(slow,az.n-1,so[0]);
	    for (iz=az.n-1; iz>0; iz--) {
		slice_get(imag,iz,qq[0]);
		LOOP( wx[imy][imx] +=
		      qq[imy][imx]; );
		
		/* upward continuation */
		slice_get(slow,iz-1,ss[0]);
		zomwex(w,iz);
		SOOP( so[ily][ilx] = ss[ily][ilx]; );
	    }

	    slice_get(imag,0,qq[0]);      /*     imaging @ iz=0 */
	    LOOP( wx[imy][imx]  += 
		  qq[imy][imx]; );		
	    
	    taper2(true,true,wx);
	    cslice_put(data,iw,wx[0]);    /* output data @ iz=0 */
	    
	} else { /* MIGRATION */
	    w = eps*aw.d - I*(aw.o+iw*aw.d); /* -1 for downward continuation */

	    cslice_get(data,iw,wx[0]);    /*  input data @ iz=0 */
	    taper2(true,true,wx);

	    slice_get(imag,0,qq[0]);      /*     imaging @ iz=0 */
	    LOOP(        qq[imy][imx] += 
		  crealf(wx[imy][imx] ); );
	    slice_put(imag,0,qq[0]);

	    slice_get(slow,0,so[0]);
	    for (iz=0; iz<az.n-1; iz++) {

		/* downward continuation */
		slice_get(slow,iz+1,ss[0]);
		zomwex(w,iz);
		SOOP( so[ily][ilx] = ss[ily][ilx]; );

		slice_get(imag,iz+1,qq[0]); /* imaging */
		LOOP(;      qq[imy][imx] += 
		     crealf(wx[imy][imx] ); );
		slice_put(imag,iz+1,qq[0]);
	    }

	} /* else */
    } /* iw */
}

/*------------------------------------------------------------*/

void cadtm(bool inv     /* forward/adjoint flag */, 
	  slice topdata /* top data [nw][nhx][nmy][nmx] */,
	  slice botdata /* bot data [nw][nhx][nmy][nmx] */)
/*< Apply upward/downward datuming >*/
{
    int iz,iw, ilx,ily;
    float complex w;

    /* loop over frequencies w */
    for (iw=0; iw<aw.n; iw++) {
	if (verb) sf_warning ("frequency %d of %d",iw+1,aw.n);

	if (inv) { /* UPWARD DATUMING */
	    w = eps*aw.d + I*(aw.o+iw*aw.d);

	    cslice_get(botdata,iw,wx[0]);
	    taper2(true,true,wx);

	    slice_get(slow,az.n-1,so[0]);
	    for (iz=az.n-1; iz>0; iz--) {
		slice_get(slow,iz-1,ss[0]);
		zomwex(w,iz);
		SOOP( so[ily][ilx] = ss[ily][ilx]; );
	    }
	    
	    taper2(true,true,wx);
	    cslice_put(topdata,iw,wx[0]);
	} else { /* DOWNWARD DATUMING */
	    w = eps*aw.d - I*(aw.o+iw*aw.d);

	    cslice_get(topdata,iw,wx[0]);
	    taper2(true,true,wx);
	    
	    slice_get(slow,0,so[0]);
	    for (iz=0; iz<az.n-1; iz++) {
		slice_get(slow,iz+1,ss[0]);
		zomwex(w,iz);
		SOOP( so[ily][ilx] = ss[ily][ilx]; );
	    }
	    
	    taper2(true,true,wx);
	    cslice_put(botdata,iw,wx[0]);
	} /* else */
    } /* iw */
}

/*------------------------------------------------------------*/

void cawfl(slice data /*      data [nw][nhx][nmy][nmx] */,
	   slice wfld /* wavefield [nw][nhx][nmy][nmx] */)
/*< Save wavefield from downward continuation >*/
{
    int iz,iw, ilx,ily;
    float complex w;

    /* loop over frequencies w */
    for (iw=0; iw<aw.n; iw++) {
	if (verb) sf_warning ("frequency %d of %d",iw+1,aw.n);
	w = eps*aw.d + I*(aw.o+iw*aw.d);

	cslice_get(data,iw,wx[0]);
	taper2(true,true,wx);
	
	taper2(true,true,wx);
	cslice_put(wfld,iw,wx[0]);

	slice_get(slow,0,so[0]);
	for (iz=0; iz<az.n-1; iz++) {	    
	    slice_get(slow,iz+1,ss[0]);
	    zomwex(w,iz);
	    SOOP( so[ily][ilx] = ss[ily][ilx]; );

	    taper2(true,true,wx);
	    cslice_put(wfld,iw,wx[0]);
	}	    
    } /* iw */
}

/*------------------------------------------------------------*/

void zomwex(float complex w,int iz)
/*< Wavefield extrapolation >*/
{
    int imy,imx,jj;
    float  s, d;
    float complex co=0, w2, cc;
    
    w2 = w*w;

    LOOP( s = so[ jy[imy] ][ jx[imx] ];
	  wx[imy][imx] *= cexpf(-w*s*az.d); );
    
    /* FFT */
    KOOP( pk[imy][imx] = 0.; );
    LOOP( pk[imy][imx] = wx[imy][imx]; );
    fft2(false,pk);

    LOOP( wx[imy][imx] = 0;
	  wt[imy][imx] = 0; );
    for (jj=0; jj<nr[iz]; jj++) {
	    
	/* w-k phase shift */
	co =       csqrtf(w2*4*sm[iz][jj]);
	KOOP( cc = csqrtf(w2*4*sm[iz][jj] + kk[imy][imx]);
	      wk[imy][imx] = 
	      pk[imy][imx] * cexpf((co-cc)*az.d); );
	
	/* IFT */
	fft2(true,wk);
	
	/* accumulate wavefield */
	LOOP( 
	    d = fabsf(4*so[ jy[imy] ][ jx[imx] ] 
		      * so[ jy[imy] ][ jx[imx] ] - 4*sm[iz][jj]);
	    d = ds2/(d*d+ds2);
	    wx[imy][imx] += wk[imy][imx]*d;
	    wt[imy][imx] += d; );
    } /* jj loop */
    LOOP( wx[imy][imx] /= wt[imy][imx]; );
    
    LOOP( s = ss[ jy[imy] ][ jx[imx] ]; 
	  wx[imy][imx] *= cexpf(-w*s*az.d); );

    taper2(true,true,wx);
}

/*------------------------------------------------------------*/

