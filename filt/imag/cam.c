/* 3-D prestack migration/modeling by split-step common-azimuth */
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

#include "cam.h"
#include "fft3.h"
#include "taper.h"
#include "slowref.h"

#include "slice.h"
/*^*/

#define LOOPxyh(a) for(ix=0;ix<ax.n;ix++){ for(iy=0;iy<ay.n;iy++){ for(ih=0;ih<ah.n;ih++){ {a} }}}
#define KOOPxyh(a) for(ix=0;ix<bx.n;ix++){ for(iy=0;iy<by.n;iy++){ for(ih=0;ih<bh.n;ih++){ {a} }}}
#define LOOPvuh(a) for(iv=0;iv<av.n;iv++){ for(iu=0;iu<au.n;iu++){ for(ih=0;ih<ah.n;ih++){ {a} }}}

#define BOUND(i,n) (i<0) ? 0 : ( (i>=n) ? n : i );
#define  KMAP(i,n) (i<n/2.) ? (i+n/2.) : (i-n/2.);

static axa az,ay,ax,aw,au,av,ah;
static axa    by,bx,         bh;
static bool verb;
static float eps;
static float  dt;
static int nrmax;

static float         ***qq; /* image */
static float complex ***pk; /* wavefield k */
static float complex ***wk; /* wavefield k */
static float complex ***wx; /* wavefield x */
static float         ***tt; /* taper */
static int           ***ms; /* multi-reference slowness map  */
static int           ***mr; /* multi-reference slowness map  */
static bool          ***skip;
static float          **ks; /* source   wavenumber */
static float          **kr; /* receiver wavenumber */
static float          **sz; /* reference slowness */
static float          **sm; /* reference slowness squared */
static float          **ss; /* slowness */
static int            **is; /* source   index */
static int            **ir; /* receiver index */
static int             *ii; /* x-line index */
static int             *ij; /* i-line index */
static int             *nr; /* number of references */

static fslice mms, mmr;

void cam_init(bool verb_,
	      float eps_,
	      float  dt_,
	      axa az_  /* depth */,
	      axa aw_  /* frequency */,
	      axa ah_  /* half-offset */,
	      axa ay_  /* i-line (data) */,
	      axa ax_  /* x-line (data) */,
	      axa au_  /* i-line (slowness/image) */,
	      axa av_  /* x-line (slowness/image) */,
	      int ntx, int nty, int nth      /* taper size */,
	      int nr1                        /* maximum number of references */,
	      int npad                       /* padding on nh */)
/*< initialize >*/
{
    int   ix, iy, ih, iu, iv;
    int       jy, jh;
    float  x,  y,  h,    k;
    float ky,     kh;

    verb=verb_;
    eps = eps_;
    dt  =  dt_;

    az = az_;
    aw = aw_;
    ay = ay_;
    ax = ax_;
    au = au_;
    av = av_;
    ah = ah_;

    bx.n = ax.n;
    bx.d =          2.0*SF_PI/(bx.n*ax.d);
    bx.o = (1==bx.n)?0:-SF_PI/      ax.d ;

    by.n = ay.n;
    by.d =          2.0*SF_PI/(by.n*ay.d);
    by.o = (1==by.n)?0:-SF_PI/      ay.d ;

    bh.n = ah.n + npad;
    bh.d =          2.0*SF_PI/(bh.n*ah.d);
    bh.d = (1==bh.n)?0:-SF_PI/      ah.d ;

    nrmax = nr1;

    fft3_init(bh.n,by.n,bx.n);

    /* allocate workspace */
    qq = sf_floatalloc3   (ah.n,au.n,av.n);  /* image */
    pk = sf_complexalloc3 (bh.n,by.n,bx.n);  /* padded wavefield */ 
    wk = sf_complexalloc3 (bh.n,by.n,bx.n);  /* k wavefield */
    wx = sf_complexalloc3 (ah.n,ay.n,ax.n);  /* x wavefield */
    tt = sf_floatalloc3   (ah.n,ay.n,ax.n);  /* taper */

    ms = sf_intalloc3     (ah.n,ay.n,ax.n);  /* MRS map */
    mr = sf_intalloc3     (ah.n,ay.n,ax.n);  /* MRS map */
    ss = sf_floatalloc2        (au.n,av.n);  /* slowness */
    sz = sf_floatalloc2       (nrmax,az.n);  /* reference slowness */
    sm = sf_floatalloc2       (nrmax,az.n);  /* reference slowness squared*/
    nr = sf_intalloc                (az.n);  /* number of reference slownesses */

    ks = sf_floatalloc2   (bh.n,ay.n);       /* source   wavenumber */
    kr = sf_floatalloc2   (bh.n,ay.n);       /* receiver wavenumber */
    is = sf_intalloc2     (ah.n,ay.n);       /* source   index */
    ir = sf_intalloc2     (ah.n,ay.n);       /* receiver index */
    ii = sf_intalloc           (ax.n);       /* midpoint index */
    ij = sf_intalloc           (ay.n);       /* midpoint index */

    skip = sf_boolalloc3 (nrmax,nrmax,az.n); /* skip S-R reference slowness combination */

    /* precompute indices */
    for (ix=0; ix<ax.n; ix++) {
	x = ax.o + ix*ax.d;

	iv = 0.5+(x-av.o)/av.d;
	iv = BOUND(iv,av.n);
	ii[ix] = iv;                     /* x-line index */
    }
    for (iy=0; iy<ay.n; iy++) {
	y = ay.o + iy*ay.d;
    
	iu = 0.5+(y-au.o)/au.d;
	iu = BOUND(iu,au.n);
	ij[iy] = iu;                     /* i-line index */

	for (ih=0; ih<ah.n; ih++) {
	    h = ah.d + ih*ah.d;
	    
	    iu = 0.5+(y-h-au.o)/au.d;
	    iu = BOUND(iu,au.n);
	    is[iy][ih] = iu;             /* source index */
	    
	    iu = 0.5+(y+h-au.o)/au.d;
	    iu = BOUND(iu,au.n);
	    ir[iy][ih] = iu;             /* receiver index */
	}
    }

    /* precompute wavenumbers */
    for (iy=0; iy<by.n; iy++) {
	jy = KMAP(iy,by.n);
	ky = by.o + jy*by.d;
	
	for (ih=0; ih<bh.n; ih++) {
	    jh = KMAP(ih,bh.n);
	    kh = bh.o + jh*bh.d;

	    k = 0.5*(ky-kh);
	    ks[iy][ih] = k*k; /* ksx^2 */

	    k = 0.5*(ky+kh);
	    kr[iy][ih] = k*k; /* krx^2 */
	}
    }    

    /* precompute taper array */
    LOOPxyh(tt[ix][iy][ih]=1.;);
    taper3(ntx,nty,nth,true,true,false,ax.n,ay.n,ah.n,tt);

    mms = fslice_init(ah.n,ax.n*ay.n,az.n,sizeof(int));
    mmr = fslice_init(ah.n,ax.n*ay.n,az.n,sizeof(int));    
}

void cam_close(void)
/*< free allocated storage >*/
{
    free(**qq); free( *qq); free( qq);
    free(**pk); free( *pk); free( pk);
    free(**wk); free( *wk); free( wk);
    free(**wx); free( *wx); free( wx);
    ;           free( *ss); free( ss);
    ;           free( *sz); free( sz);
    ;           free( *sm); free( sm);
    ;           ;           free( nr);
    ;           free( *ks); free( ks);
    ;           free( *kr); free( kr);
    ;           free( *is); free( is);
    ;           free( *ir); free( ir);
    free(**tt); free( *tt); free( tt);
    free(**ms); free( *ms); free( ms);
    free(**mr); free( *mr); free( mr);
    
    free(**skip); free( *skip); free(skip);

    fslice_close(mms);
    fslice_close(mmr);
    
    fft3_close();
}

void cam( bool inv   /* migration/modeling flag */, 
	  slice data /* data       [nw][nx][ny][nh] */,
	  slice imag /* image file [nz][nv][nu][nh] */,
	  slice slow /* slowness   [nz][nv][nu]     */)
/*< Apply migration/modeling >*/
{
    int iz,iw,ix,iy,ih,iu,iv, j,js,jr;
    int       jx;
    float sy, kx;
    float complex cshift=0, cref=0, w, w2, cs, cr, kh, kss, krr;

    if (!inv) { /* prepare image for migration */
	LOOPvuh( qq[iv][iu][ih] = 0.0; );
	for (iz=0; iz<az.n; iz++) {
	    slice_put(imag,iz,qq[0][0]);
	}
    }

    /* compute reference slowness */
    for (iz=0; iz<az.n; iz++) {
	slice_get(slow,iz,ss[0]);
	nr[iz] = slowref(nrmax,dt/az.d,au.n*av.n,ss[0],sz[iz],sm[iz]);
	if (verb) sf_warning("nr[%d]=%d",iz,nr[iz]);

	/* create MRS map */
	LOOPxyh( ms[ix][iy][ih] = 0; 
		 mr[ix][iy][ih] = 0.; );
	for (j=0; j<nr[iz]-1; j++) {
	    sy = 0.5*(sz[iz][j]+sz[iz][j+1]);
	    LOOPxyh( if(ss[ ii[ix] ][ is[iy][ih] ] > sy) ms[ix][iy][ih]++;
		     if(ss[ ii[ix] ][ ir[iy][ih] ] > sy) mr[ix][iy][ih]++; );
	}
	fslice_put(mms,iz,ms[0][0]);
	fslice_put(mmr,iz,mr[0][0]);
	for (js=0; js<nr[iz]; js++) {
	    for (jr=0; jr<nr[iz]; jr++) {
		skip[iz][js][jr] = true;
	    }
	}
	LOOPxyh( skip[iz][ ms[ix][iy][ih] ][ mr[ix][iy][ih] ] = false; );
    }
    for (iz=0; iz<az.n-1; iz++) {
	for (j=0; j<nr[iz]; j++) {
	    sm[iz][j] = 0.5*(sm[iz][j]+sm[iz+1][j]);
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
	    LOOPxyh( wx[ix][iy][ih] = qq[ ii[ix] ][ ij[iy] ][ih]; );
	  
	    /* loop over migrated depths z */
	    for (iz=az.n-2; iz>=0; iz--) {

		/* w-x @ bottom */
		LOOPxyh( sy = 0.5*(ss[ ii[ix] ][ is[iy][ih] ] + 
				   ss[ ii[ix] ][ ir[iy][ih] ]);
			 cshift = cexpf(-w*sy*az.d);
			 wx[ix][iy][ih] *= cshift; );
		LOOPxyh( pk[ix][iy][ih] *= tt[ix][iy][ih] ; );

		/* FFT */
		KOOPxyh( pk[ix][iy][ih] = 0.;);
		LOOPxyh( pk[ix][iy][ih] = wx[ix][iy][ih]; );
		fft3(false,pk);		
		
		fslice_get(mms,iz,ms[0][0]);
		fslice_get(mmr,iz,mr[0][0]);

		LOOPxyh( wx[ix][iy][ih] = 0; );
		for (js=0; js<nr[iz]; js++) {
		    for (jr=0; jr<nr[iz]; jr++) {
			if (skip[iz][js][jr]) continue; /* skip S-R reference combinations */

			/* w-k phase shift */
			cref = csqrtf(w2*sm[iz][js]) + csqrtf(w2*sm[iz][jr]); /* w so */
			KOOPxyh( jx = KMAP(ix,bx.n);
				 kx = bx.o + jx*bx.d; 
				 cs = csqrtf(w2*sm[iz][js]+ks[iy][ih]);
				 cr = csqrtf(w2*sm[iz][jr]+kr[iy][ih]);
				 kh = kx*(cr-cs)/(cr+cs); /* comaz approximation */
				 kss = 0.5*(kx-kh);
				 krr = 0.5*(kx+kh);
				 kss = kss*kss + ks[iy][ih];
				 krr = krr*krr + kr[iy][ih];
				 cs = csqrtf(w2*sm[iz][js] + kss);
				 cr = csqrtf(w2*sm[iz][jr] + krr);
				 cshift = cexpf((cref-cs-cr)*az.d);  /* w so - kzs - kzr */
				 wk[ix][iy][ih] = 
				 pk[ix][iy][ih] * cshift; ); 
			
			/* IFT */
			fft3(true,wk);

			/* accumulate wavefield */
			LOOPxyh( if (ms[ix][iy][ih]==js && 
				     mr[ix][iy][ih]==jr ) 
				 wx[ix][iy][ih] += wk[ix][iy][ih]; );
		    } /* jr loop */
		} /* js loop */

		/* w-x at top */
		slice_get(slow,iz,ss[0]);
		LOOPxyh( sy = 0.5*(ss[ ii[ix] ][ is[iy][ih] ] + 
				   ss[ ii[ix] ][ ir[iy][ih] ]);
			 cshift = cexpf(-w*sy*az.d);
			 wx[ix][iy][ih] *= cshift; ); 

		/* imaging condition */
		slice_get(imag,iz,qq[0][0]);
		LOOPxyh( wx   [ix]    [iy] [ih] +=
			 qq[ii[ix]][ij[iy]][ih]; );

	    } /* iz */
	    
	    LOOPxyh( wx[ix][iy][ih] *= tt[ix][iy][ih]; ); /* taper wavefield */
	    cslice_put(data,iw,wx[0][0]);    /* put wavefield = data on file */
	    
	} else { /* MIGRATION */
	    slice_get(slow,0,ss[0]);                      /* slowness at top */

	    cslice_get(data,iw,wx[0][0]);    /* get wavefield = data on file */
	    LOOPxyh( wx[ix][iy][ih] *= tt[ix][iy][ih]; ); /* taper wavefield */

	    /* loop over migrated depths z */
	    for (iz=0; iz< az.n-1; iz++) {

		/* imaging condition */
		slice_get(imag,iz,qq[0][0]);
		LOOPxyh(        qq[ii[ix]][ij[iy]][ih] += 
			 crealf(wx   [ix]    [iy] [ih] ); );
		slice_put(imag,iz,qq[0][0]);

		/* w-x @ top */
		LOOPxyh( sy = 0.5*(ss[ ii[ix] ][ is[iy][ih] ] + 
				   ss[ ii[ix] ][ ir[iy][ih] ]);
			 cshift = conjf(cexpf(-w*sy*az.d));
			 wx[ix][iy][ih] *= cshift; );

		/* FFT */
		KOOPxyh( pk[ix][iy][ih] = 0.;);
		LOOPxyh( pk[ix][iy][ih] = wx[ix][iy][ih]; );
		fft3(false,pk);

		fslice_get(mms,iz,ms[0][0]);
		fslice_get(mmr,iz,mr[0][0]);
		
		LOOPxyh( wx[ix][iy][ih] = 0; );
		for (js=0; js<nr[iz]; js++) {
		    for (jr=0; jr<nr[iz]; jr++) {
			if (skip[iz][js][jr]) continue;
		
			/* w-k phase shift */
			cref= csqrtf(w2*sm[iz][js]) 
			    + csqrtf(w2*sm[iz][jr]); /* w so */
			KOOPxyh( jx = KMAP(ix,bx.n);
				 kx = bx.o + jx*bx.d; 
				 cs = csqrtf(w2*sm[iz][js] + ks[iy][ih]);
				 cr = csqrtf(w2*sm[iz][jr] + kr[iy][ih]);
				 kh = kx*(cr-cs)/(cr+cs); /* comaz approximation */
				 kss = 0.5*(kx-kh);
				 krr = 0.5*(kx+kh);
				 kss = kss*kss + ks[iy][ih];
				 krr = krr*krr + kr[iy][ih];
				 cs = csqrtf(w2*sm[iz][js] + kss);
				 cr = csqrtf(w2*sm[iz][jr] + krr);
				 cshift = conjf(cexpf((cref-cs-cr)*az.d)); /* w so - kzs - kzr */
				 wk[ix][iy][ih] = 
				 pk[ix][iy][ih] * cshift; );
 
			/* IFT */
			fft3(true,wk);

			/* accumulate wavefield */
			LOOPxyh( if (ms[ix][iy][ih]==js && 
				     mr[ix][iy][ih]==jr ) 
				 wx[ix][iy][ih] += wk[ix][iy][ih]; );
		    } /* jr loop */
		} /* js loop */

		/* w-x @ bottom */
		slice_get(slow,iz+1,ss[0]);
		LOOPxyh( sy = 0.5*(ss[ ii[ix] ][ is[iy][ih] ] + 
				   ss[ ii[ix] ][ ir[iy][ih] ]);
			 cshift = conjf(cexpf(-w*sy*az.d));
			 wx[ix][iy][ih] *= cshift; );
		LOOPxyh( wx[ix][iy][ih] *= tt[ix][iy][ih]; );
	    } /* iz */
	    
	    /* imaging condition @ bottom */
	    slice_get(imag,az.n-1,qq[0][0]);
	    LOOPxyh(        qq[ii[ix]][ij[iy]][ih] += 
		     crealf(wx   [ix]    [iy] [ih]); );
	    slice_put(imag,az.n-1,qq[0][0]);
	} /* else */
    } /* iw */
}

