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
#include "taper.h"
#include "slowref.h"

#include "slice.h"
/*^*/

#define LOOPxy(a) for(ix=0;ix<nx;ix++){ for(iy=0;iy<ny;iy++){ {a} }}

static int nx,ny,nz,nrmax;
static float dz;
static float         **ss; /* slowness */
static float         **sz; /* reference slowness */
static float         **sm; /* reference slowness squared */
static int            *nr; /* number of references */

static float         **qq; /* image */
static float         **kk; /* wavenumber */
static float         **tt; /* taper */

static float complex **pp; /* wavefield */
static float complex **wt; /* wavefield top */
static float complex **wb; /* wavefield bot */

static int           **mm; /* multi-reference slowness map  */
static float         **ma; /* multi-reference slowness mask */

void split2_init(int nz1, float dz1  /* depth */,
		 int ny1, float dy1  /* in-line midpoint */,
		 int nx1, float dx1  /* cross-line midpoint */,
		 int ntx, int nty    /* taper size */,
		 int nr1             /* maximum number of references */) 
/*< initialize >*/
{
    int ix, iy, jx, jy;
    float x0, y0, dx, dy, kx, ky;

    nz = nz1;
    dz = dz1;

    nx = nx1;
    dx =        2.0*SF_PI/(nx*dx1);
    x0 = (1==nx)?0:-SF_PI/dx1;

    ny = ny1;
    dy =        2.0*SF_PI/(ny*dy1);
    y0 = (1==ny)?0:-SF_PI/dy1;

    nrmax = nr1;

    fft2_init(ny,nx);

    /* allocate workspace */
    ss = sf_floatalloc2 (ny,nx);    /* slowness */
    sz = sf_floatalloc2 (nrmax,nz); /* reference slowness */
    sm = sf_floatalloc2 (nrmax,nz); /* reference slowness squared*/
    nr = sf_intalloc    (      nz); /* number of reference slownesses */

    qq = sf_floatalloc2   (ny,nx);  /* image */
    kk = sf_floatalloc2   (ny,nx);  /* wavenumber */
    tt = sf_floatalloc2   (ny,nx);  /* taper */

    pp = sf_complexalloc2 (ny,nx);  /* wavefield */
    wt = sf_complexalloc2 (ny,nx);  /* wavefield top */
    wb = sf_complexalloc2 (ny,nx);  /* wavefield bot */

    mm = sf_intalloc2     (ny,nx);  /* MRS map */
    ma = sf_floatalloc2   (ny,nx);  /* MRS mask */

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

    /* precompute taper array */
    LOOPxy(tt[ix][iy]=1.;);
    taper2(ntx,nty,true,true,nx,ny,tt);
}

void split2_close(void)
/*< free allocated storage >*/
{
    fft2_close();

    free(*pp); free(pp);
    free(*wt); free(wt);
    free(*wb); free(wb);

    free(*ss); free(ss);
    free(*sz); free(sz);
    free(*sm); free(sm);
    free( nr);

    free(*qq); free(qq);
    free(*kk); free(kk);
    free(*tt); free(tt);

    free(*mm); free(mm);
    free(*ma); free(ma);
}

void split2(bool verb                   /* verbosity flag */, 
	    bool inv                    /* migration/modeling flag */, 
	    float eps                   /* stability factor */,  
	    int nw, float dw, float w0  /* frequency (radian) */,
	    sf_file data                /* data  [nw][nx][ny] */,
	    slice imag                  /* image  [nz][nx][ny] */,
	    slice slow                  /* slowness  [nz][nx][ny] */,
	    float dt                    /* time error */)
/*< Apply migration/modeling >*/
{
    int iz,iw,ix,iy,ir,nxy;
    float complex cshift, cref=0, w, w2;

    if (!inv) { /* prepare image for migration */
	LOOPxy( qq[ix][iy] = 0.0; );
	for (iz=0; iz<nz; iz++) {
	    slice_put(imag,iz,qq[0]);
	}
    }

    /* compute reference slowness */
    nxy = nx*ny;
    for (iz=0; iz<nz; iz++) {
	slice_get(slow,iz,  ss[0]);
	nr[iz] = slowref(nrmax,dt/dz,nxy,ss[0],sz[iz],sm[iz]);
	if (verb) sf_warning("nr[%d]=%d",iz,nr[iz]);
    }
    for (iz=0; iz<nz-1; iz++) {
	for (ir=0; ir<nr[iz]; ir++) {
	    sm[iz][ir] = 2*(sm[iz][ir]+sm[iz+1][ir]);
	}
    }

    /* loop over frequencies w */
    for (iw=0; iw<nw; iw++) {
	if (verb) sf_warning ("frequency %d of %d",iw+1,nw);

	w = eps*dw + I*(w0+iw*dw);
	w2 = w*w;

	if (inv) { /* MODELING */
	    /* start from bottom */
	    slice_get(slow,nz-1,ss[0]);
	    
	    /* imaging condition */
	    slice_get(imag,nz-1,qq[0]);
	    LOOPxy( pp[ix][iy] = qq[ix][iy]; );

	    /* loop over migrated depths z */
	    for (iz=nz-2; iz>=0; iz--) {
		
		/* w-x @ bottom */
		LOOPxy( cshift = cexpf(-w*ss[ix][iy]*dz)*tt[ix][iy];
			pp[ix][iy] *= cshift; );
		
		/* FFT */
		fft2(false,pp);

		slice_get(slow,iz,ss[0]);

		/* create MRS map */
		LOOPxy( mm[ix][iy] = 0; );
		for (ir=0; ir<nr[iz]-1; ir++) {
		    LOOPxy( if(ss[ix][iy] > 0.5*(sz[iz][ir]+sz[iz][ir+1])) mm[ix][iy]++; );
		}

		LOOPxy( wt[ix][iy] = 0; );
		for (ir=0; ir<nr[iz]; ir++) {
		    /* w-k phase shift */
		    cref = csqrtf(w2*sm[iz][ir]);
		    LOOPxy( cshift = cexpf((cref-csqrtf(w2*sm[iz][ir]+kk[ix][iy]))*dz);
			    wb[ix][iy] = pp[ix][iy]*cshift; );

		    /* IFT */
		    fft2(true,wb);

		    /* create MRS mask */
		    LOOPxy( ma[ix][iy]= (mm[ix][iy]==ir)?1.:0.; );

		    /* accumulate wavefield */
		    LOOPxy( wt[ix][iy] += wb[ix][iy] * ma[ix][iy]; );
		    
		} /* ir loop */

		slice_get(imag,iz,qq[0]);
		
		/* w-x at top */
		LOOPxy( cshift = cexpf(-w*ss[ix][iy]*dz);
			pp[ix][iy] = qq[ix][iy] + wt[ix][iy]*cshift; );
	    } /* iz */

	    /* taper */
	    LOOPxy( pp[ix][iy] *= tt[ix][iy]; );

	    sf_complexwrite(pp[0],nxy,data);

	} else { /* MIGRATION */

	    sf_complexread(pp[0],nxy,data);

	    /* taper */
	    LOOPxy( pp[ix][iy] *= tt[ix][iy]; );

	    slice_get(slow,0,ss[0]);


	    /* loop over migrated depths z */
	    for (iz=0; iz< nz-1; iz++) {

		/* imaging condition */
		slice_get(imag,iz,qq[0]);
		LOOPxy( qq[ix][iy] += crealf(pp[ix][iy]); );
		slice_put(imag,iz,qq[0]);

		/* w-x @ top */
		LOOPxy( cshift= cexpf(-w*ss[ix][iy]*dz);
			pp[ix][iy] *= conjf(cshift); );

		/* FFT */
		fft2(false,pp);
		
		slice_get(slow,iz+1,ss[0]);

		/* create MRS map */
		LOOPxy( mm[ix][iy] = 0; );
		for (ir=0; ir<nr[iz]-1; ir++) {
		    LOOPxy( if(ss[ix][iy] > 0.5*(sz[iz+1][ir]+sz[iz+1][ir+1])) 
			    mm[ix][iy]++; );
		}

		LOOPxy( wb[ix][iy] = 0; );
		for (ir=0; ir<nr[iz]; ir++) {
		    /* w-k phase shift */
		    cref = csqrtf(w2*sm[iz][ir]);
		    LOOPxy( cshift = cexpf((cref-csqrtf(w2*sm[iz][ir]+kk[ix][iy]))*dz);
			    wt[ix][iy] = pp[ix][iy]*conjf(cshift); );
		    
		    /* IFT */
		    fft2(true,wt);

		    /* create MRS mask */
		    LOOPxy( ma[ix][iy]=(mm[ix][iy]==ir)?1.:0.; );

		    /* accumulate wavefield */
		    LOOPxy( wb[ix][iy] += wt[ix][iy] * ma[ix][iy]; );
		    
		} /* ir loop */
		
		/* w-x @ bottom */
		LOOPxy( cshift = cexpf(-w*ss[ix][iy]*dz)*tt[ix][iy];
			pp[ix][iy] = wb[ix][iy]*conjf(cshift); );
	    } /* iz */

	    /* imaging condition @ bottom */ 
	    slice_get(imag,nz-1,qq[0]);
	    LOOPxy( qq[ix][iy] += crealf(pp[ix][iy]); );
	    slice_put(imag,nz-1,qq[0]);

	} /* else */
    } /* iw */
}

