/* 2-D prestack migration/modeling by split-step DSR */
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

#include "dsr2.h"
#include "fft2.h"
#include "taper.h"
#include "slowref.h"

#define LOOPxh(a)  for(ix=0;ix<nx;ix++){ for(ih=0;ih<nh; ih++){ {a} }}
#define LOOPxh2(a) for(ix=0;ix<nx;ix++){ for(ih=0;ih<nh2;ih++){ {a} }}
#define LOOPuh(a)  for(iu=0;iu<nu;iu++){ for(ih=0;ih<nh; ih++){ {a} }}

#include "slice.h"
/*^*/

static int nx,nh,nh2,nz,nu,nrmax;
static float dz;

static float         **qq;             /* image */
static int     **is, **ir, *ii;        /* indices */
static float   **ks, **kr;             /* wavenumber */

static float         **sz;             /* reference slowness */
static float         **sm;             /* reference slowness squared */
static int            *nr;             /* number of references */

static float complex **pk;   /* wavefield */
static float complex **wk;   /* wavefield k */
static float complex **wx;   /* wavefield x */

static int           **ms, **mr; /* multi-reference slowness map  */
static bool          ***skip;
static float         **ma; /* multi-reference slowness mask */
static fslice mms, mmr;

void dsr2_init(int nz1, float dz1             /* depth */,
	       int nh1, float dh1, float h01  /* half-offset */,
	       int nx1, float dx1, float x01  /* midpoint */,
	       int nu1, float du,  float u0   /* slowness grid */,
	       int ntx, int nth               /* taper size */,
	       int nr1                        /* number of references */,
	       int npad                       /* padding on nh */)
/*< initialize >*/
{
    int   ix,  ih,  iu;
    int   jx,  jh;
    float  x,   h;
    float  kx, kh, k;
    float  dx, dh;
    float   x0, h0;
    
    nz = nz1;
    dz = dz1;

    nx = nx1;
    dx = 2.0*SF_PI/(nx*dx1);
    x0 =    -SF_PI/    dx1 ;

    nh = nh1;
    nh2 = nh+npad;

    dh = 2.0*SF_PI/(nh2*dh1);
    h0 =    -SF_PI/     dh1 ;

    nu = nu1;

    nrmax = nr1;

    fft2_init(nh2,nx);

    /* allocate workspace */
    sz = sf_floatalloc2  (nrmax,nz);       /* reference slowness */
    sm = sf_floatalloc2  (nrmax,nz);       /* reference slowness squared*/
    nr = sf_intalloc     (      nz);       /* number of reference slownesses */

    qq = sf_floatalloc2     (nh,nu);       /* image */

    ks = sf_floatalloc2   (nh2,nx);        /* source wavenumber */
    kr = sf_floatalloc2   (nh2,nx);        /* receiver wavenumber */
    is = sf_intalloc2     (nh, nx);        /* source reference */
    ir = sf_intalloc2     (nh, nx);        /* receiver reference */
    ii = sf_intalloc          (nx);        /* midpoint reference */

    pk = sf_complexalloc2 (nh2,nx);        /* padded wavefield */ 
    wx = sf_complexalloc2 (nh, nx);        /* x wavefield */
    wk = sf_complexalloc2 (nh2,nx);        /* k wavefield */

    ms = sf_intalloc2     (nh, nx);        /* MRS map source */
    mr = sf_intalloc2     (nh, nx);        /* MRS map receiver */
    ma = sf_floatalloc2   (nh, nx);        /* MRS mask */

    skip = sf_boolalloc3 (nrmax,nrmax,nz); /* skip slowness combination */

    /* precompute wavenumbers */
    for (ix=0; ix<nx; ix++) {
	jx = (ix < nx/2)? ix + nx/2: ix - nx/2;
	kx = x0 + jx*dx;
	x = x01 + ix*dx1;
	
	iu = 0.5+(x-u0)/du;
	if      (iu <   0) iu=0;
	else if (iu >= nu) iu=nu-1;
	ii[ix] = iu;

	for (ih=0; ih<nh; ih++) {
	    h = h01 + ih*dh1;
	    
	    iu = 0.5+(x-h-u0)/du;
	    if      (iu <   0) iu=0;
	    else if (iu >= nu) iu=nu-1;
	    is[ix][ih] = iu;

	    iu = 0.5+(x+h-u0)/du;
	    if      (iu <   0) iu=0;
	    else if (iu >= nu) iu=nu-1;
	    ir[ix][ih] = iu;
	}

	for (ih=0; ih<nh2; ih++) {
	    jh = (ih < nh2/2)? ih + nh2/2: ih - nh2/2;
	    kh = h0 + jh*dh;

	    k = 0.5*(kx-kh);
	    ks[ix][ih] = k*k;

	    k = 0.5*(kx+kh);
	    kr[ix][ih] = k*k;
	}
    }

    /* precompute taper array */
    taper2_init(ntx,nth);

    mms = fslice_init(nh,nx,nz,sizeof(int));
    mmr = fslice_init(nh,nx,nz,sizeof(int));
}

void dsr2_close(void)
/*< free allocated storage >*/
{
    free( *pk); free( pk);
    free( *wk); free( wk);
    free( *wx); free( wx);

    free( *sz); free( sz);
    free( *sm); free( sm);
    free(  nr);

    free( *qq); free( qq);
    free( *ks); free( ks);
    free( *kr); free( kr);
    free( *is); free( is);
    free( *ir); free( ir);
    free(  ii);

    free( *ms); free( ms);
    free( *mr); free( mr);
    free( *ma); free( ma);
    
    free(**skip); 
    free( *skip); free(skip);
    fslice_close(mms);
    fslice_close(mmr);
    taper2_close();
    fft2_close();
}

void dsr2(bool verb                   /* verbosity flag */, 
	  bool inv                    /* migration/modeling flag */, 
	  float eps                   /* stability factor */,  
	  int nw, float dw, float w0  /* frequency (radian) */,
	  sf_file data                /* data       [nw][nx][nh] */,
	  slice imag                  /* image file [nz][nu][nh] */,
	  float **slow                /* slowness   [nz][nx]     */,
	  float dt                    /* time error */)
/*< Apply migration/modeling >*/
{
    int iz,iw,ix,ih,iu, j,k;
    float sy, *si;
    float complex cshift, w, w2, cs, cr, cref;

    if (!inv) { /* prepare image for migration */
	LOOPuh( qq[iu][ih] = 0.0; );
	for (iz=0; iz<nz; iz++) {
	    slice_put(imag,iz,qq[0]);
	}
    }

    /* compute reference slowness */
    for (iz=0; iz<nz; iz++) {
	si = slow[iz];
	nr[iz] = slowref(nrmax,dt/dz,nu,si,sz[iz],sm[iz]);
	if (verb) sf_warning("nr[%d]=%d",iz,nr[iz]);
	/* create MRS map */
	LOOPxh( ms[ix][ih] = 0; mr[ix][ih] = 0.;);
	for (j=0; j<nr[iz]-1; j++) {
	    sy = 0.5*(sz[iz][j]+sz[iz][j+1]);
	    LOOPxh( if(si[ is[ix][ih] ] > sy) ms[ix][ih]++;
		    if(si[ ir[ix][ih] ] > sy) mr[ix][ih]++; );
	}
	fslice_put(mms,iz,ms[0]);
	fslice_put(mmr,iz,mr[0]);
	for (j=0; j<nr[iz]; j++) {
	    for (k=0; k<nr[iz]; k++) {
		skip[iz][j][k] = true;
	    }
	}
	LOOPxh( skip[iz][ms[ix][ih]][mr[ix][ih]] = false; );
    }
    for (iz=0; iz<nz-1; iz++) {
	for (j=0; j<nr[iz]; j++) {
	    sm[iz][j] = 0.5*(sm[iz][j]+sm[iz+1][j]);
	}
    }
    
    /* loop over frequencies w */
    for (iw=0; iw<nw; iw++) {
	if (verb) sf_warning ("frequency %d of %d",iw+1,nw);

	w = eps*dw + I*(w0+iw*dw);
	w2 = w*w;

	if (inv) { /* MODELING */
	    /* start from bottom */
	    si = slow[nz-1];

	    /* imaging condition */
	    slice_get(imag,nz-1,qq[0]);	    
	    LOOPxh( wx[ix][ih] = qq[ ii[ix] ][ih];  );

	    /* loop over migrated depths z */
	    for (iz=nz-2; iz>=0; iz--) {
		/* w-x @ bottom */
		LOOPxh2( pk[ix][ih] = 0.;);
		taper2(true,false,nx,nh,wx);
		LOOPxh( sy = 0.5*(si[ is[ix][ih] ] + 
				  si[ ir[ix][ih] ]);
			cshift = cexpf(-w*sy*dz);
			pk[ix][ih] = 
			wx[ix][ih] * cshift; );
		
		/* FFT */
		fft2(false,pk);		
		
		si = slow[iz];
		fslice_get(mms,iz,ms[0]);
		fslice_get(mmr,iz,mr[0]);

		LOOPxh( wx[ix][ih] = 0; );
		for (j=0; j<nr[iz]; j++) {
		    for (k=0; k<nr[iz]; k++) {
			if (skip[iz][j][k]) continue; /* skip S-R reference combinations */

			/* w-k phase shift */
			cref =       csqrtf(w2*sm[iz][j])
			    +        csqrtf(w2*sm[iz][k]);
			LOOPxh2(cs = csqrtf(w2*sm[iz][j] + ks[ix][ih]);
				cr = csqrtf(w2*sm[iz][k] + kr[ix][ih]);
				cshift = cexpf((cref-cs-cr)*dz); 
				wk[ix][ih] = 
				pk[ix][ih]*cshift; ); 
		
			/* IFT */
			fft2(true,wk);

			/* create MRS mask */
			LOOPxh( ma[ix][ih]= (ms[ix][ih]==j && 
					     mr[ix][ih]==k)?1.:0.; );

			/* accumulate wavefield */
			LOOPxh( wx[ix][ih] += wk[ix][ih] * ma[ix][ih]; );
		    }
		}
		
		slice_get(imag,iz,qq[0]);

		/* w-x at top */
		LOOPxh( sy = 0.5*(si[ is[ix][ih] ] + 
				  si[ ir[ix][ih] ]);
			cshift = cexpf(-w*sy*dz);
			wx   [ix] [ih] = 
			qq[ii[ix]][ih] + 
			wx   [ix] [ih] * cshift; );
	    } /* iz */

	    taper2(true,false,nx,nh,wx);
	    sf_complexwrite(wx[0],nx*nh,data);

	} else { /* MIGRATION */
	    si = slow[0];

	    sf_complexread(wx[0],nx*nh,data);
	    taper2(true,false,nx,nh,wx);

	    /* loop over migrated depths z */
	    for (iz=0; iz< nz-1; iz++) {

		/* imaging condition */
		slice_get(imag,iz,qq[0]);
		LOOPxh(        qq[ii[ix]][ih] += 
			crealf(wx   [ix] [ih] ); );
		slice_put(imag,iz,qq[0]);

		/* w-x @ top */
		LOOPxh2( pk[ix][ih] = 0.;);
		LOOPxh( sy = 0.5*(si[ is[ix][ih] ] + 
				  si[ ir[ix][ih] ]);
			cshift = conjf(cexpf(-w*sy*dz));
			pk[ix][ih] = 
			wx[ix][ih] * cshift; );

		/* FFT */
		fft2(false,pk);

		si = slow[iz+1];
		fslice_get(mms,iz,ms[0]);
		fslice_get(mmr,iz,mr[0]);
		
		LOOPxh( wx[ix][ih] = 0; );
		for (j=0; j<nr[iz]; j++) {
		    for (k=0; k<nr[iz]; k++) {
			if (skip[iz][j][k]) continue;
		
			/* w-k phase shift */
			cref = csqrtf(w2*sm[iz][j])+csqrtf(w2*sm[iz][k]);
			LOOPxh2( cs = csqrtf(w2*sm[iz][j] + ks[ix][ih]);
				 cr = csqrtf(w2*sm[iz][k] + kr[ix][ih]);
				 cshift = conjf(cexpf((cref-cs-cr)*dz)); 
				 wk[ix][ih] = 
				 pk[ix][ih] * cshift; ); 
		    
			
			/* IFT */
			fft2(true,wk);

			/* create MRS mask */
			LOOPxh( ma[ix][ih]= (ms[ix][ih]==j && 
					     mr[ix][ih]==k)?1.:0.; );

			/* accumulate wavefield */
			LOOPxh( wx[ix][ih] += wk[ix][ih] * ma[ix][ih]; );
		    }
		} /* j loop */


		/* w-x @ bottom */
		LOOPxh( sy = 0.5*(si[ is[ix][ih] ] + 
				  si[ ir[ix][ih] ]);
			cshift = conjf(cexpf(-w*sy*dz));
			wx[ix][ih] *= cshift; );
		taper2(true,false,nx,nh,wx);
	    } /* iz */
	    
	    /* imaging condition @ bottom */
	    slice_get(imag,nz-1,qq[0]);
	    LOOPxh(        qq[ii[ix]][ih] += 
		    crealf(wx   [ix] [ih]); );
	    slice_put(imag,nz-1,qq[0]);
	} /* else */
    } /* iw */
}
