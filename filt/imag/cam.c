/* 3-D common-azimuth migration/modeling using split-step */
/*q
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

#include "cam.h"
#include "fft3.h"
#include "taper.h"
#include "slowref.h"

#include "slice.h"
/*^*/

#define LOOP(a) for(ihx=0;ihx<ahx.n;ihx++){ for(imy=0;imy<amy.n;imy++){ for(imx=0;imx<amx.n;imx++){ {a} }}}
#define KOOP(a) for(ihx=0;ihx<bhx.n;ihx++){ for(imy=0;imy<bmy.n;imy++){ for(imx=0;imx<bmx.n;imx++){ {a} }}}

#define INDEX(x,a) 0.5+(x-a.o)/a.d;
#define BOUND(i,n) (i<0) ? 0 : ( (i>=n) ? n : i );
#define  KMAP(i,n) (i<n/2.) ? (i+n/2.) : (i-n/2.);
#define X2K(a,b,p) b.n=a.n+p; b.d=2.0*SF_PI/(b.n*a.d); b.o=(1==b.n)?0:-SF_PI/a.d;

static axa az,amx,amy,aw,alx,aly,ahx;
static axa    bmx,bmy,           bhx;
static bool verb;
static float eps;
static float ds, ds2;

static float         ***qq; /* image */
static float complex ***pk; /* wavefield k */
static float complex ***wk; /* wavefield k */
static float complex ***wx; /* wavefield x */
static float         ***wt; /* interpolation weight */
static float          **ksx;/* source   wavenumber  */
static float          **krx;/* receiver wavenumber  */
static float          **sz; /* reference slowness   */
static float          **sm; /* reference slowness squared */
static float          **ss; /* slowness */
static int            **is; /* source   index */
static int            **ir; /* receiver index */
static int             *jx; /* i-line index */
static int             *jy; /* x-line index */
static int             *nr; /* number of references */

void cam_init(bool verb_,
	      float eps_,
	      float  dt,
	      axa az_   /* depth */,
	      axa aw_   /* frequency */,
	      axa ahx_  /* half-offset */,
	      axa amx_  /* i-line (data) */,
	      axa amy_  /* x-line (data) */,
	      axa alx_  /* i-line (slowness/image) */,
	      axa aly_  /* x-line (slowness/image) */,
	      int ntx, int nty, int nth      /* taper size */,
	      int nrmax                      /* maximum number of references */,
	      int padmx,int padmy,int padhx  /* padding in the k domain */,
	      slice slow)
/*< initialize >*/
{
    int   imy, imx, ihx, ilx, ily, iz, j;
    int        jmx, jhx;
    float  my,  mx,  hx,    k;
    float      kmx, khx;

    verb=verb_;
    eps = eps_;

    az = az_;
    aw = aw_;
    amx= amx_;
    amy= amy_;
    alx= alx_;
    aly= aly_;
    ahx= ahx_;

    ds  = dt/az.d;
    ds2 = ds*ds;
    ds2 *= ds2;

    /* construct K-domain axes */
    X2K(amx,bmx,padmx);
    X2K(amy,bmy,padmy);
    X2K(ahx,bhx,padhx);
    fft3_init(bmx.n,bmy.n,bhx.n);

    /* allocate workspace */
    ss = sf_floatalloc2   (alx.n,aly.n      );  /* slowness */
    sz = sf_floatalloc2          (nrmax,az.n);  /* reference slowness */
    sm = sf_floatalloc2          (nrmax,az.n);  /* reference slowness squared*/
    nr = sf_intalloc                   (az.n);  /* number of reference slownesses */

    qq = sf_floatalloc3   (amx.n,amy.n,ahx.n);  /* image */
    wx = sf_complexalloc3 (amx.n,amy.n,ahx.n);  /* x wavefield */
    wt = sf_floatalloc3   (amx.n,amy.n,ahx.n);  /* interpolation weight */

    wk = sf_complexalloc3 (bmx.n,bmy.n,bhx.n);  /* k wavefield */
    pk = sf_complexalloc3 (bmx.n,bmy.n,bhx.n);  /* padded wavefield */ 

    ksx= sf_floatalloc2   (bmx.n,      bhx.n);  /* source   wavenumber */
    krx= sf_floatalloc2   (bmx.n,      bhx.n);  /* receiver wavenumber */

    jy = sf_intalloc      (      amy.n      );  /* midpoint index */
    jx = sf_intalloc      (amx.n            );  /* midpoint index */
    is = sf_intalloc2     (amx.n,      ahx.n);  /* source   index */
    ir = sf_intalloc2     (amx.n,      ahx.n);  /* receiver index */

    /* precompute indices */
    for (imy=0; imy<amy.n; imy++) {
	my = amy.o + imy*amy.d;
	ily     = INDEX( my,aly);
	jy[imy] = BOUND(ily,aly.n);                     /* x-line index */
    }
    for (imx=0; imx<amx.n; imx++) {
	mx = amx.o + imx*amx.d;
	ilx     = INDEX( mx,alx);
	jx[imx] = BOUND(ilx,alx.n);                     /* i-line index */

	for (ihx=0; ihx<ahx.n; ihx++) {
	    hx = ahx.o + ihx*ahx.d;
	    
	    ilx          = INDEX(mx-hx,alx);
	    is[ihx][imx] = BOUND(  ilx,alx.n);          /* source index */

	    ilx          = INDEX(mx+hx,alx);
	    ir[ihx][imx] = BOUND(  ilx,alx.n);          /* receiver index */
	}
    }

    /* precompute wavenumbers */
    for (imx=0; imx<bmx.n; imx++) {
	jmx = KMAP(imx,bmx.n);
	kmx = bmx.o + jmx*bmx.d;
	
	for (ihx=0; ihx<bhx.n; ihx++) {
	    jhx = KMAP(ihx,bhx.n);
	    khx = bhx.o + jhx*bhx.d;

	    k = 0.5*(kmx-khx);
	    ksx[ihx][imx] = k*k; /* ksx^2 */

	    k = 0.5*(kmx+khx);
	    krx[ihx][imx] = k*k; /* krx^2 */
	}
    }    

    /* precompute taper array */
    taper3_init(nth,nty,ntx);

    /* compute reference slowness */
    for (iz=0; iz<az.n; iz++) {
	slice_get(slow,iz,ss[0]);
	nr[iz] = slowref(nrmax,ds,alx.n*aly.n,ss[0],sz[iz],sm[iz]);
	if (verb) sf_warning("nr[%d]=%d",iz,nr[iz]);
    }
    for (iz=0; iz<az.n-1; iz++) {
	for (j=0; j<nr[iz]; j++) {
	    sm[iz][j] = 0.5*(sm[iz][j]+sm[iz+1][j]);
	}
    }
}

void cam_close(void)
/*< free allocated storage >*/
{
    free(**qq); free( *qq); free( qq);
    free(**pk); free( *pk); free( pk);
    free(**wk); free( *wk); free( wk);
    free(**wx); free( *wx); free( wx);
    free(**wt); free( *wt); free( wt);
    ;           free( *ss); free( ss);
    ;           free( *sz); free( sz);
    ;           free( *sm); free( sm);
    ;           ;           free( nr);
    ;           free( *ksx);free( ksx);
    ;           free( *krx);free( krx);
    ;           free( *is); free( is);
    ;           free( *ir); free( ir);
    
    fft3_close();
    taper3_close();
}

void cam( bool inv   /* migration/modeling flag */, 
	  slice data /* data       [nw][nhx][nmy][nmx] */,
	  slice imag /* image file [nz][nhx][nly][nlx] */,
	  slice slow /* slowness   [nz]     [nly][nlx]     */)
/*< Apply migration/modeling >*/
{
    int iz,iw,imy,imx,ihx,js,jr;
    int       jmy;
    float sy, kmy;
    float complex cshift=0, cref=0, w, w2, cs, cr, khy, kss, krr;
    float d, dsc,drc;

    if (!inv) { /* prepare image for migration */
	LOOP( qq[ihx][imy][imx] = 0.0; );
	for (iz=0; iz<az.n; iz++) {
	    slice_put(imag,iz,qq[0][0]);
	}
    }
    
    /* loop over frequencies w */
    for (iw=0; iw<aw.n; iw++) {
	if (verb) sf_warning ("frequency %d of %d",iw+1,aw.n);

	w = eps*aw.d + I*(aw.o+iw*aw.d);
	w2 = w*w;

	if (inv) { /* MODELING */
	    slice_get(slow,az.n-1,ss[0]);              /* slowness at bottom */

	    /* imaging condition @ bottom */
	    slice_get(imag,az.n-1,qq[0][0]);
	    LOOP( wx[ihx][imy][imx] = qq[ihx][imy][imx]; );
	  
	    /* loop over migrated depths z */
	    for (iz=az.n-2; iz>=0; iz--) {

		/* w-x @ bottom */
		LOOP( sy = 0.5*(ss[ jy[imy] ][ is[ihx][imx] ] + 
				   ss[ jy[imy] ][ ir[ihx][imx] ]);
			 cshift = cexpf(-w*sy*az.d);
			 wx[ihx][imy][imx] *= cshift; );
		taper3(true,true,false,ahx.n,amy.n,amx.n,wx);

		/* FFT */
		KOOP( pk[ihx][imy][imx] = 0.;);
		LOOP( pk[ihx][imy][imx] = wx[ihx][imy][imx]; );
		fft3(false,pk);		
		
		LOOP( wx[ihx][imy][imx] = 0; );
		LOOP( wt[ihx][imy][imx] = 0; );
		for (js=0; js<nr[iz]; js++) {
		    for (jr=0; jr<nr[iz]; jr++) {

			/* w-k phase shift */
			cref = csqrtf(w2*sm[iz][js]) + csqrtf(w2*sm[iz][jr]); /* w so */
			KOOP( jmy = KMAP(imy,bmy.n);
			      kmy = bmy.o + jmy*bmy.d; 
			      cs  = csqrtf(w2*sm[iz][js]+ksx[ihx][imx]);
			      cr  = csqrtf(w2*sm[iz][jr]+krx[ihx][imx]);
			      khy = kmy*(cr-cs)/(cr+cs); /* comaz approximation */
			      kss = 0.5*(kmy-khy);
			      krr = 0.5*(kmy+khy);
			      kss = kss*kss + ksx[ihx][imx];
			      krr = krr*krr + krx[ihx][imx];
			      cs  = csqrtf(w2*sm[iz][js] + kss);
			      cr  = csqrtf(w2*sm[iz][jr] + krr);
			      cshift = cexpf((cref-cs-cr)*az.d);  /* w so - kzs - kzr */
			      wk[ihx][imy][imx] = 
			      pk[ihx][imy][imx] * cshift; ); 
			
			/* IFT */
			fft3(true,wk);

			/* accumulate wavefield */
			LOOP( 
			    dsc = fabsf( ss[ jy[imy] ][ is[ihx][imx] ] *
					 ss[ jy[imy] ][ is[ihx][imx] ] - sm[iz][js]);
			    drc = fabsf( ss[ jy[imy] ][ ir[ihx][imx] ] *
					 ss[ jy[imy] ][ ir[ihx][imx] ] - sm[iz][jr]);
			    d = sqrt(dsc*dsc + drc*drc);
			    d = 1./(d*d+ds2);
			    wx[ihx][imy][imx] += wk[ihx][imy][imx]*d;
			    wt[ihx][imy][imx] += d; );			
		    } /* jr loop */
		} /* js loop */
		LOOP( wx[ihx][imy][imx] /= wt[ihx][imy][imx]; );
		
		/* w-x at top */
		slice_get(slow,iz,ss[0]);
		LOOP( sy = 0.5*(ss[ jy[imy] ][ is[ihx][imx] ] + 
				   ss[ jy[imy] ][ ir[ihx][imx] ]);
			 cshift = cexpf(-w*sy*az.d);
			 wx[ihx][imy][imx] *= cshift; ); 

		/* imaging condition */
		slice_get(imag,iz,qq[0][0]);
		LOOP( wx[ihx][imy][imx]  +=
			 qq[ihx][imy][imx]; );

	    } /* iz */

	    taper3(true,true,false,ahx.n,amy.n,amx.n,wx);
	    cslice_put(data,iw,wx[0][0]);    /* put wavefield = data on file */
	    
	} else { /* MIGRATION */
	    slice_get(slow,0,ss[0]);                      /* slowness at top */

	    cslice_get(data,iw,wx[0][0]);    /* get wavefield = data on file */
	    taper3(true,true,false,ahx.n,amy.n,amx.n,wx);

	    /* loop over migrated depths z */
	    for (iz=0; iz< az.n-1; iz++) {

		/* imaging condition */
		slice_get(imag,iz,qq[0][0]);
		LOOP(       qq[ihx][imy][imx] += 
		     crealf(wx[ihx][imy][imx] ); );
		slice_put(imag,iz,qq[0][0]);

		/* w-x @ top */
		LOOP( sy = 0.5*(ss[ jy[imy] ][ is[ihx][imx] ] + 
				ss[ jy[imy] ][ ir[ihx][imx] ]);
			 cshift = conjf(cexpf(-w*sy*az.d));
			 wx[ihx][imy][imx] *= cshift; );

		/* FFT */
		KOOP( pk[ihx][imy][imx] = 0.; );
		LOOP( pk[ihx][imy][imx] = wx[ihx][imy][imx]; );
		fft3(false,pk);

		LOOP( wx[ihx][imy][imx] = 0; );
		LOOP( wt[ihx][imy][imx] = 0; );
		for (js=0; js<nr[iz]; js++) {
		    for (jr=0; jr<nr[iz]; jr++) {

			/* w-k phase shift */
			cref= csqrtf(w2*sm[iz][js]) 
			    + csqrtf(w2*sm[iz][jr]); /* w so */
			KOOP( jmy = KMAP(imy,bmy.n);
			      kmy = bmy.o + jmy*bmy.d; 
			      cs  = csqrtf(w2*sm[iz][js] + ksx[ihx][imx]);
			      cr  = csqrtf(w2*sm[iz][jr] + krx[ihx][imx]);
			      khy = kmy*(cr-cs)/(cr+cs); /* comaz approximation */
			      kss = 0.5*(kmy-khy);
			      krr = 0.5*(kmy+khy);
			      kss = kss*kss + ksx[ihx][imx];
			      krr = krr*krr + krx[ihx][imx];
			      cs  = csqrtf(w2*sm[iz][js] + kss);
			      cr  = csqrtf(w2*sm[iz][jr] + krr);
			      cshift = conjf(cexpf((cref-cs-cr)*az.d)); /* w so - kzs - kzr */
			      wk[ihx][imy][imx] = 
			      pk[ihx][imy][imx] * cshift; );
			
			/* IFT */
			fft3(true,wk);

			/* accumulate wavefield */
			LOOP( 
			    dsc = fabsf( ss[ jy[imy] ][ is[ihx][imx] ] *
					 ss[ jy[imy] ][ is[ihx][imx] ] - sm[iz][js]);
			    drc = fabsf( ss[ jy[imy] ][ ir[ihx][imx] ] *
					 ss[ jy[imy] ][ ir[ihx][imx] ] - sm[iz][jr]);
			    d = sqrt(dsc*dsc + drc*drc);
			    d = 1./(d*d+ds2);
			    wx[ihx][imy][imx] += wk[ihx][imy][imx]*d;
			    wt[ihx][imy][imx] += d; );
		    } /* jr loop */
		} /* js loop */
		LOOP( wx[ihx][imy][imx] /= wt[ihx][imy][imx]; );

		/* w-x @ bottom */
		slice_get(slow,iz+1,ss[0]);
		LOOP( sy = 0.5*(ss[ jy[imy] ][ is[ihx][imx] ] + 
				ss[ jy[imy] ][ ir[ihx][imx] ]);
			 cshift = conjf(cexpf(-w*sy*az.d));
			 wx[ihx][imy][imx] *= cshift; );
		taper3(true,true,false,ahx.n,amy.n,amx.n,wx);
	    } /* iz */
	    
	    /* imaging condition @ bottom */
	    slice_get(imag,az.n-1,qq[0][0]);
	    LOOP(qq[ihx][imy][imx] +=  crealf(wx[ihx][imy][imx]); );
	    slice_put(imag,az.n-1,qq[0][0]);
	} /* else */
    } /* iw */
}

