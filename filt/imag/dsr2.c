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

#include "slice.h"
/*^*/

static int nx, nh, nz, ny, **is, **ir, **ii;
static float dz, **qq, **ks, **kr, *s;
static float complex **pp;

void dsr2_init(int nz1, float dz1            /* depth */,
	       int nh1, float dh1, float h01 /* half-offset */,
	       int nx1, float dx1, float x01 /* midpoint */,
	       int ny1, float dy, float y0   /* slowness grid */)
/*< initialize >*/
{
    int ix, ih, jx, jh, iy;
    float x, h, x0, h0, dx, dh, kx, kh, k;

    nz = nz1;
    dz = dz1;

    nx = nx1;
    dx = 2.0*SF_PI/(nx*dx1);
    x0 = -SF_PI/dx1;

    nh = nh1;
    dh = 2.0*SF_PI/(nh*dh1);
    h0 = -SF_PI/dh1;

    ny = ny1;

    /* allocate workspace */
    pp  = sf_complexalloc2 (nh,nx);
    fft2_init(nh,nx);

    s = sf_floatalloc (nz);      /* reference slowness */
    qq = sf_floatalloc2 (nh,ny); /* image */
    ks = sf_floatalloc2 (nh,nx); /* source wavenumber */
    kr = sf_floatalloc2 (nh,nx); /* receiver wavenumber */
    is = sf_intalloc2(nh,nx);    /* source reference */
    ir = sf_intalloc2(nh,nx);    /* receiver reference */
    ii = sf_intalloc2(nh,nx);    /* midpoint reference */

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

	    iy = 0.5+(x-h-y0)/dy;
	    if (iy < 0) iy=0;
	    else if (iy >= ny) iy=ny-1;
	    is[ix][ih] = iy;

	    iy = 0.5+(x+h-y0)/dy;
	    if (iy < 0) iy=0;
	    else if (iy >= ny) iy=ny-1;
	    ir[ix][ih] = iy;

	    iy = 0.5+(x-y0)/dy;
	    if (iy < 0) iy=0;
	    else if (iy >= ny) iy=ny-1;
	    ii[ix][ih] = iy;
	}
    }    
}

void dsr2_close(void)
/*< free allocated storage >*/
{
    free(*pp);
    free(pp);
    free(s);
    free(*qq);
    free(qq);
    free(*ks);
    free(ks);
    free(*kr);
    free(kr);
    free(*is);
    free(is);
    free(*ir);
    free(ir);
}

void dsr2(bool verb                   /* verbosity flag */, 
	  bool inv                    /* migration/modeling flag */, 
	  float eps                   /* stability factor */,  
	  int nw, float dw, float w0  /* frequency (radian) */,
	  float complex *** cp        /* data [nw][nx][nh] */,
	  slice imag                  /* image file [nz][ny][nh] */,
	  float **slow                /* slowness [nz][nx] */)
/*< Apply migration/modeling >*/
{
    int iz,iw,ix,ih,iy;
    float sy, *si;
    float complex cshift, cref, w, w2, **pp, cs, cr;

    if (!inv) { /* prepare image for migration */
	for (iy=0; iy<ny; iy++) {      
	    for (ih=0; ih<nh; ih++) {
		qq[iy][ih] = 0.0;
	    }
	}
	for (iz=0; iz<nz; iz++) {
	    slice_put(imag,iz,qq[0]);
	}
    }

    /* compute reference slowness squared */
    for (iz=0; iz<nz; iz++) {
	/* median */
	s[iz] = sf_quantile(0.5*ny,ny,slow[iz]);
	s[iz] *= s[iz];
    }	
    /* take average in layer */
    for (iz=0; iz<nz-1; iz++) {
	s[iz] = 0.5*(s[iz]+s[iz+1]);
    }

    /* loop over frequencies w */
    for (iw=0; iw<nw; iw++) {
	if (verb) sf_warning ("frequency %d of %d",iw+1,nw);

	w = eps*dw + I*(w0+iw*dw);
	w2 = w*w;

	pp = cp[iw];

	if (inv) { /* modeling */
	    /* start from bottom */
	    slice_get(imag,nz-1,qq[0]);
	    si = slow[nz-1];
	    
	    for (ix=0; ix<nx; ix++) {
		for (ih=0; ih<nh; ih++) {
		    pp[ix][ih] = qq[ii[ix][ih]][ih];
		}
	    }

	    /* loop over migrated depths z */
	    for (iz=nz-2; iz>=0; iz--) {

		/* space-domain, part 1 */
		for (ix=0; ix<nx; ix++) {
		    for (ih=0; ih<nh; ih++) {
			sy = si[is[ix][ih]]+si[ir[ix][ih]];
			cshift = cexpf(-0.5*w*sy*dz);
			pp[ix][ih] *= cshift; /* add tapering later */
		    }
		}

		fft2(false,pp);		

		/* phase shift */
		cref = 2.*csqrtf(w2*s[iz]);
		for (ix=0; ix<nx; ix++) {
		    for (ih=0; ih<nh; ih++) {
			cs = csqrtf(w2*s[iz]+ks[ix][ih]);
			cr = csqrtf(w2*s[iz]+kr[ix][ih]);
			cshift = cexpf((cref-cs-cr)*dz); 
			pp[ix][ih] *= cshift; 
		    }
		}
		
		fft2(true,pp);
		
		slice_get(imag,iz,qq[0]);
		si = slow[iz];
		
		/* space-domain, part 1 */
		for (ix=0; ix<nx; ix++) {
		    for (ih=0; ih<nh; ih++) {
			sy = si[is[ix][ih]]+si[ir[ix][ih]];
			cshift = cexpf(-0.5*w*sy*dz);
			pp[ix][ih] = qq[ii[ix][ih]][ih] + pp[ix][ih]*cshift; 
                        /* add tapering later */
		    }
		}
	    } /* iz */
	} else {
	    /* loop over migrated depths z */
	    si = slow[0];

	    for (iz=0; iz< nz-1; iz++) {
		slice_get(imag,iz,qq[0]);

		for (ix=0; ix<nx; ix++) {
		    for (ih=0; ih<nh; ih++) {
			qq[ii[ix][ih]][ih] += crealf(pp[ix][ih]); /* imaging cond. */
			sy = si[is[ix][ih]]+si[ir[ix][ih]];
			cshift = conjf(cexpf(-0.5*w*sy*dz));
			pp[ix][ih] *= cshift;
		    }
		}
		slice_put(imag,iz,qq[0]);

		fft2(false,pp);
		
		/* phase shift */
		cref = csqrtf(w2*s[iz]);
		for (ix=0; ix<nx; ix++) {
		    for (ih=0; ih<nh; ih++) {
			cs = csqrtf(w2*s[iz]+ks[ix][ih]);
			cr = csqrtf(w2*s[iz]+kr[ix][ih]);
			cshift = conjf(cexpf((cref-cs-cr)*dz)); 
			pp[ix][ih] *= cshift; 
		    }
		}
		
		fft2(true,pp);

		si = slow[iz+1];

		for (ix=0; ix<nx; ix++) {
		    for (ih=0; ih<nh; ih++) {
			sy = si[is[ix][ih]]+si[ir[ix][ih]];
			cshift = conjf(cexpf(-0.5*w*sy*dz));
			pp[ix][ih] *= cshift;
		    }
		}
	    } /* iz */
	    
	    /* arrive to bottom */
	    slice_get(imag,nz-1,qq[0]);
	    
	    for (ix=0; ix<nx; ix++) {
		for (ih=0; ih<nh; ih++) {
		    qq[ii[ix][ih]][ih] += crealf(pp[ix][ih]); /* imaging condition */ 
		}
	    }	    
	    slice_put(imag,nz-1,qq[0]);
	} /* else */
    } /* iw */
}
