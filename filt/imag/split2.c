/* 3-D split-step migration/modeling */
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

#include "split2.h"
#include "fft2.h"

#include "slice.h"
/*^*/

static int nx, ny, nz;
static float dz, **ss, **qq, **kk, *s;
static float complex **pp;

void split2_init(int nz1, float dz1 /* depth */,
		 int ny1, float dy1  /* in-line midpoint */,
		 int nx1, float dx1  /* cross-line midpoint */)
/*< initialize >*/
{
    int ix, iy, jx, jy;
    float x0, y0, dx, dy, kx, ky;

    nz = nz1;
    dz = dz1;

    nx = nx1;
    dx = 2.0*SF_PI/(nx*dx1);
    x0 = -SF_PI/dx1;

    ny = ny1;
    dy = 2.0*SF_PI/(ny*dy1);
    y0 = -SF_PI/dy1;

    /* allocate workspace */
    pp  = sf_complexalloc2 (ny,nx);
    fft2_init(ny,nx);

    ss = sf_floatalloc2 (ny,nx); /* slowness */
    s = sf_floatalloc (nz); /* reference slowness */
    qq = sf_floatalloc2 (ny,nx); /* image */
    kk = sf_floatalloc2 (ny,nx); /* wavenumber */

    /* precompute wavenumbers */
    for (ix=0; ix<nx; ix++) {
	jx = (ix < nx/2)? ix + nx/2: ix - nx/2;
	kx = x0 + jx*dx;
	kx *= kx;
	
	for (iy=0; iy<ny; iy++) {
	    jy = (iy < ny/2)? iy + ny/2: iy - ny/2;
	    ky = y0 + jy*dy;
	    kk[ix][iy] = ky*ky + kx;
	}
    }    
}

void split2_close(void)
/*< free allocated storage >*/
{
    free(*pp);
    free(pp);
    free(*ss);
    free(ss);
    free(s);
    free(*qq);
    free(qq);
    free(*kk);
    free(kk);
}

void split2(bool verb                   /* verbosity flag */, 
	    bool inv                    /* migration/modeling flag */, 
	    float eps                   /* stability factor */,  
	    int nw, float dw, float w0  /* frequency (radian) */,
	    float complex *** cp        /* data [nw][nx][ny] */,
	    slice imag                  /* image [nz][nx][ny] */,
	    slice slow                  /* file with slowness */)
/*< Apply migration/modeling >*/
{
    int iz,iw,ix,iy;
    float complex cshift, cref, w, w2, **pp;

    if (!inv) { /* prepare image for migration */
	for (ix=0; ix<nx; ix++) {      
	    for (iy=0; iy<ny; iy++) {
		qq[ix][iy] = 0.0;
	    }
	}
	for (iz=0; iz<nz; iz++) {
	    slice_put(imag,iz,qq[0]);
	}
    }

    
    for (iz=0; iz<nz; iz++) {
	/* compute reference slowness squared */
	slice_get(slow,iz,ss[0]);
	/* median */
	s[iz] = sf_quantile(0.5*nx*ny,nx*ny,ss[0]);
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
	    slice_get(slow,nz-1,ss[0]);
	    
	    for (ix=0; ix<nx; ix++) {
		for (iy=0; iy<ny; iy++) {
		    pp[ix][iy] = qq[ix][iy];
		}
	    }

	    /* loop over migrated depths z */
	    for (iz=nz-2; iz>=0; iz--) {

		/* space-domain, part 1 */
		for (ix=0; ix<nx; ix++) {
		    for (iy=0; iy<ny; iy++) {
			cshift = cexpf(-0.5*w*ss[ix][iy]*dz);
			pp[ix][iy] *= cshift; /* add tapering later */
		    }
		}

		fft2(false,pp);		

		/* phase shift */
		cref = csqrtf(w2*s[iz]);
		for (ix=0; ix<nx; ix++) {
		    for (iy=0; iy<ny; iy++) {
			cshift = cexpf((cref-csqrtf(w2*s[iz]+kk[ix][iy]))*dz); 
			pp[ix][iy] *= cshift; 
		    }
		}

		
		fft2(true,pp);

		
		slice_get(imag,iz,qq[0]);
		slice_get(slow,iz,ss[0]);
		
		/* space-domain, part 1 */
		for (ix=0; ix<nx; ix++) {
		    for (iy=0; iy<ny; iy++) {
			cshift = cexpf(-0.5*w*ss[ix][iy]*dz);
			pp[ix][iy] = qq[ix][iy] + pp[ix][iy]*cshift; /* add tapering later */
		    }
		}
	    } /* iz */
	} else {
	    /* loop over migrated depths z */
	    for (iz=0; iz< nz-1; iz++) {
		slice_get(imag,iz,qq[0]);
		slice_get(slow,iz,ss[0]);


		for (ix=0; ix<nx; ix++) {
		    for (iy=0; iy<ny; iy++) {
			qq[ix][iy] += crealf(pp[ix][iy]); /* imaging condition */
			cshift = conjf(cexpf(-0.5*w*ss[ix][iy]*dz));
			pp[ix][iy] *= cshift;
		    }
		}
		slice_put(imag,iz,qq[0]);
		

		fft2(false,pp);

		/* phase shift */
		cref = csqrtf(w2*s[iz]);
		for (ix=0; ix<nx; ix++) {
		    for (iy=0; iy<ny; iy++) {
			cshift = conjf(cexpf((cref-csqrtf(w2*s[iz]+kk[ix][iy]))*dz)); 
			pp[ix][iy] *= cshift; 
		    }
		}

		
		fft2(true,pp);

		for (ix=0; ix<nx; ix++) {
		    for (iy=0; iy<ny; iy++) {
			cshift = conjf(cexpf(-0.5*w*ss[ix][iy]*dz));
			pp[ix][iy] *= cshift;
		    }
		}
	    } /* iz */

	    
	    /* arrive to bottom */
	    slice_get(imag,nz-1,qq[0]);
	    slice_get(slow,nz-1,ss[0]);
	    
	    for (ix=0; ix<nx; ix++) {
		for (iy=0; iy<ny; iy++) {
		    qq[ix][iy] += crealf(pp[ix][iy]); /* imaging condition */ 
		}
	    }	    
	    slice_put(imag,nz-1,qq[0]);
	} /* else */
    } /* iw */
}
