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

#include "comaz.h"
#include "fft2.h"
#include "taper.h"
#include "slowref.h"

#define LOOPxh(a) for(ix=0;ix<nx;ix++){ for(ih=0;ih<nh;ih++){ {a} }}
#define LOOPuh(a) for(iu=0;iu<nu;iu++){ for(ih=0;ih<nh;ih++){ {a} }}

#include "slice.h"
/*^*/

static int nx,nh,nz,nu,nrmax;
static float dz;

static float         **qq;             /* image */
static int     **is, **ir, **ii;       /* indices */
static float   **ks, **kr;             /* wavenumber */
static float         **tt;             /* taper */

static float         **sz;             /* reference slowness */
static float         **sm;             /* reference slowness squared */
static int            *nr;             /* number of references */

static float complex **pp; /* wavefield */
static float complex **wt; /* wavefield top */
static float complex **wb; /* wavefield bot */

static int           **mm; /* multi-reference slowness map  */
static float         **ma; /* multi-reference slowness mask */

void comaz_init(int nz1, float dz1             /* depth */,
	       int nh1, float dh1, float h01  /* half-offset */,
	       int nx1, float dx1, float x01  /* midpoint */,
	       int nu1, float du,  float u0   /* slowness grid */,
	       int nt                         /* taper size */,
	       int nr1                        /* maximum number of references */)
/*< initialize >*/
{
    int ix, ih, jx, jh, iu;
    float x, h, x0, h0, dx, dh, kx, kh, k;

    nz = nz1;
    dz = dz1;

    nx = nx1;
    dx = 2.0*SF_PI/(nx*dx1);
    x0 =    -SF_PI/dx1;

    nh = nh1;
    dh = 2.0*SF_PI/(nh*dh1);
    h0 =    -SF_PI/dh1;

    nu = nu1;

    nrmax = nr1;

    fft2_init(nh,nx);

    /* allocate workspace */
    sz = sf_floatalloc2  (nrmax,nz); /* reference slowness */
    sm = sf_floatalloc2  (nrmax,nz); /* reference slowness squared*/
    nr = sf_intalloc     (      nz); /* number of reference slownesses */

    qq = sf_floatalloc2   (nh,nu);   /* image */

    ks = sf_floatalloc2   (nh,nx);   /* source wavenumber */
    kr = sf_floatalloc2   (nh,nx);   /* receiver wavenumber */
    is = sf_intalloc2     (nh,nx);   /* source reference */
    ir = sf_intalloc2     (nh,nx);   /* receiver reference */
    ii = sf_intalloc2     (nh,nx);   /* midpoint reference */

    tt = sf_floatalloc2   (nh,nx);   /* taper */

    pp = sf_complexalloc2 (nh,nx);   /* wavefield */ 
    wt = sf_complexalloc2 (nh,nx);  /* wavefield top */
    wb = sf_complexalloc2 (nh,nx);  /* wavefield bot */

    mm = sf_intalloc2     (nh,nx);  /* MRS map */
    ma = sf_floatalloc2   (nh,nx);  /* MRS mask */

    /* precompute wavenumbers */
    for (ix=0; ix<nx; ix++) {
	jx = (ix < nx/2)? ix + nx/2: ix - nx/2;
	kx = x0 + jx*dx;
	x = x01 + ix*dx1;

	for (ih=0; ih<nh; ih++) {
	    jh = (ih < nh/2)? ih + nh/2: ih - nh/2;
	    kh = h0 + jh*dh;
	    h = h01 + ih*dh1;

	    k = 0.5*(kx-kh);
	    ks[ix][ih] = k*k;

	    k = 0.5*(kx+kh);
	    kr[ix][ih] = k*k;

	    iu = 0.5+(x-h-u0)/du;
	    if      (iu <   0) iu=0;
	    else if (iu >= nu) iu=nu-1;
	    is[ix][ih] = iu;

	    iu = 0.5+(x+h-u0)/du;
	    if      (iu <   0) iu=0;
	    else if (iu >= nu) iu=nu-1;
	    ir[ix][ih] = iu;

	    iu = 0.5+(x-u0)/du;
	    if      (iu <   0) iu=0;
	    else if (iu >= nu) iu=nu-1;
	    ii[ix][ih] = iu;
	}
    }    

    /* precompute taper array */
    LOOPxh(tt[ix][ih]=1.;);
/*    taper2(nt,0,nx,nh,tt); */
}

void comaz_close(void)
/*< free allocated storage >*/
{
    free(*pp); free(pp);
    free(*wt); free(wt);
    free(*wb); free(wb);

    free(*sz); free(sz);
    free(*sm); free(sm);
    free( nr);

    free(*qq); free(qq);
    free(*ks); free(ks);
    free(*kr); free(kr);
    free(*is); free(is);
    free(*ir); free(ir);
    
    free(*tt); free(tt);

    free(*mm); free(mm);
    free(*ma); free(ma);
    
}

void comaz(bool verb                   /* verbosity flag */, 
	   bool inv                    /* migration/modeling flag */, 
	   float eps                   /* stability factor */,  
	   int nw, float dw, float w0  /* frequency (radian) */,
	   sf_file data                /* data       [nw][nx][nh] */,
	   slice imag                  /* image file [nz][nu][nh] */,
	   float **slow                /* slowness   [nz][nx]     */,
	   float dt                    /* time error */)
/*< Apply migration/modeling >*/
{
    int iz,iw,ix,ih,iu, jr;
    float sy, *si;
    float complex cshift, cref, w, w2, cs, cr;

    if (!inv) { /* prepare image for migration */
	LOOPuh( qq[iu][ih] = 0.0; );
	for (iz=0; iz<nz; iz++) {
	    slice_put(imag,iz,qq[0]);
	}
    }

    /* compute reference slowness */
    for (iz=0; iz<nz; iz++) {
	nr[iz] = slowref(nrmax,dt/dz,nu,slow[iz],sz[iz],sm[iz]);
	if (verb) sf_warning("nr[%d]=%d",iz,nr[iz]);
    }
    for (iz=0; iz<nz-1; iz++) {
	for (jr=0; jr<nr[iz]; jr++) {
	    sm[iz][jr] = 0.5*(sm[iz][jr]+sm[iz+1][jr]);
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
	    LOOPxh( pp[ix][ih] = qq[ii[ix][ih]][ih];  );

	    /* loop over migrated depths z */
	    for (iz=nz-2; iz>=0; iz--) {
		/* w-x @ bottom */
		LOOPxh( sy = 0.5*(si[ is[ix][ih] ] + si[ ir[ix][ih] ]);
			cshift = cexpf(-w*sy*dz);
			pp[ix][ih] *= cshift; );
		
		/* FFT */
		fft2(false,pp);		
		
		si = slow[iz];
		
		/* create MRS map */
		LOOPxh( mm[ix][ih] = 0; );

		LOOPxh( wt[ix][ih] = 0; );

		for (jr=0; jr<nr[iz]; jr++) {
		    /* w-k phase shift */
		    cref = 2.*csqrtf(   w2*sm[iz][jr]);
		    LOOPxh(	cs = csqrtf(w2*sm[iz][jr]+ks[ix][ih]);
				cr = csqrtf(w2*sm[iz][jr]+kr[ix][ih]);
				cshift = cexpf((cref-cs-cr)*dz); 
				wb[ix][ih] = pp[ix][ih]*cshift; ); 
		
		    /* IFT */
		    fft2(true,wb);

		    /* create MRS mask */
		    LOOPxh( ma[ix][ih]= (mm[ix][ih]==ir)?1.:0.; );

		    /* accumulate wavefield */
		    LOOPxh( wt[ix][ih] += wb[ix][ih] * ma[ix][ih]; );
		}

		
		slice_get(imag,iz,qq[0]);

		/* w-x at top */
		LOOPxh( sy = 0.5*(si[ is[ix][ih] ] + si[ ir[ix][ih] ]);
			cshift = cexpf(-w*sy*dz);
			pp[ix][ih] = qq[ii[ix][ih]][ih] + wt[ix][ih]*cshift; );
	    } /* iz */
	    
	    /* taper */
	    LOOPxh( pp[ix][ih] *= tt[ix][ih]; );

	    sf_complexwrite(pp[0],nx*nh,data);

	} else { /* MIGRATION */
	    si = slow[0];

	    sf_complexread(pp[0],nx*nh,data);

	    /* taper */
	    LOOPxh( pp[ix][ih] *= tt[ix][ih]; );

	    /* loop over migrated depths z */
	    for (iz=0; iz< nz-1; iz++) {

		/* imaging condition */
		slice_get(imag,iz,qq[0]);
		LOOPxh( qq[ii[ix][ih]][ih] += crealf(pp[ix][ih]); );
		slice_put(imag,iz,qq[0]);

		/* w-x @ top */
		LOOPxh( sy = 0.5*(si[ is[ix][ih] ] + si[ ir[ix][ih] ]);
			cshift = conjf(cexpf(-w*sy*dz));
			pp[ix][ih] *= cshift; );

		/* FFT */
		fft2(false,pp);

		si = slow[iz+1];
		
		/* create MRS map */
		LOOPxh( mm[ix][ih] = 0; );
		for (jr=0; jr<nr[iz]-1; jr++) {
		    LOOPxh( if(si[ ii[ix][ih] ] > 0.5*(sz[iz+1][jr]+sz[iz+1][jr+1])) 
			    mm[ix][ih]++; );
		}

		LOOPxh( wb[ix][ih] = 0; );
		for (jr=0; jr<nr[iz]; jr++) {
		
		    /* w-k phase shift */
		    cref = 2.*csqrtf(   w2*sm[iz][jr]);
		    LOOPxh( cs = csqrtf(w2*sm[iz][jr]+ks[ix][ih]);
			    cr = csqrtf(w2*sm[iz][jr]+kr[ix][ih]);
			    cshift = conjf(cexpf((cref-cs-cr)*dz)); 
			    wt[ix][ih] = pp[ix][ih]*cshift; ); 
		    

		    /* IFT */
		    fft2(true,wt);

		    /* create MRS mask */
		    LOOPxh( ma[ix][ih]= (mm[ix][ih]==jr)?1.:0.; );

		    /* accumulate wavefield */
		    LOOPxh( wb[ix][ih] += wt[ix][ih] * ma[ix][ih]; );
		    
		} /* jr loop */


		/* w-x @ bottom */
		LOOPxh( sy = 0.5*(si[ is[ix][ih] ] + si[ ir[ix][ih] ]);
			cshift = conjf(cexpf(-w*sy*dz));
			pp[ix][ih] = wb[ix][ih]*cshift; );
	    } /* iz */
	    
	    /* imaging condition @ bottom */
	    slice_get(imag,nz-1,qq[0]);
	    LOOPxh( qq[ii[ix][ih]][ih] += crealf(pp[ix][ih]); );
	    slice_put(imag,nz-1,qq[0]);
	} /* else */
    } /* iw */
}
