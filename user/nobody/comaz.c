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
#include "fft3.h"
#include "taper.h"
#include "slowref.h"

#define LOOPxyh(a)  for(ix=0;ix<nx;ix++){ for(iy=0;iy<ny;iy++){ for(ih=0;ih<nh; ih++){ {a} }}}
#define LOOPxyh2(a) for(ix=0;ix<nx;ix++){ for(iy=0;iy<ny;iy++){ for(ih=0;ih<nh2;ih++){ {a} }}}
#define LOOPvuh(a)  for(iv=0;iv<nv;iv++){ for(iu=0;iu<nu;iu++){ for(ih=0;ih<nh; ih++){ {a} }}}

#include "slice.h"
/*^*/

static int nx,ny,nh,nh2,nz,nu,nv,nrmax;
static float dz,x0,dx;

static float         ***qq;            /* image */
static float   **ss;                   /* slowness */
static int     **is, **ir, *ii, *ij;   /* indices */
static float   **ks, **kr;             /* wavenumber */

static float         **sz;             /* reference slowness */
static float         **sm;             /* reference slowness squared */
static int            *nr;             /* number of references */

static float complex ***pk; /* wavefield */
static float complex ***wk; /* wavefield k */
static float complex ***wx; /* wavefield x */

static int    ***ms, ***mr; /* multi-reference slowness map  */
static bool   ***skip;
static fslice mms, mmr;

void comaz_init(int nz1, float dz1             /* depth */,
		int nh1, float dh1, float h01  /* half-offset */,
		int ny1, float dy1, float y01  /* i-line (data) */,
		int nx1, float dx1, float x01  /* x-line (data) */,
		int nu1, float du,  float u0   /* i-line (slowness/image) */,
		int nv1, float dv,  float v0   /* x-line (slowness/image) */,
		int ntx, int nty, int nth      /* taper size */,
		int nr1                        /* number of references */,
		int npad                       /* padding on nh */)
/*< initialize >*/
{
    int ix, iy, jy, ih, jh, iu, iv;
    float x, y, h, y0, h0, dy, dh, ky, kh, k;

    nz = nz1;
    dz = dz1;

    nx = nx1;
    dx = 2.0*SF_PI/(nx*dx1);
    x0 =    -SF_PI/    dx1 ;

    ny = ny1;
    dy = 2.0*SF_PI/(ny*dy1);
    y0 =    -SF_PI/    dy1 ;

    nh = nh1;
    nh2 = nh+npad;

    dh = 2.0*SF_PI/(nh2*dh1);
    h0 =    -SF_PI/     dh1 ;

    nu = nu1;
    nv = nv1;

    nrmax = nr1;

    fft3_init(nh2,ny,nx);

    /* allocate workspace */
    sz = sf_floatalloc2  (nrmax,nz);      /* reference slowness */
    sm = sf_floatalloc2  (nrmax,nz);      /* reference slowness squared*/
    nr = sf_intalloc     (      nz);      /* number of reference slownesses */

    qq = sf_floatalloc3  (nh,nu,nv);      /* image */
    ss = sf_floatalloc2     (nu,nv);      /* slowness */

    ks = sf_floatalloc2   (nh2,ny);       /* source wavenumber */
    kr = sf_floatalloc2   (nh2,ny);       /* receiver wavenumber */
    is = sf_intalloc2     (nh ,ny);       /* source reference */
    ir = sf_intalloc2     (nh ,ny);       /* receiver reference */
    ii = sf_intalloc          (nx);       /* midpoint reference */
    ij = sf_intalloc          (ny);       /* midpoint reference */

    pk = sf_complexalloc3 (nh2,ny, nx);   /* padded wavefield */ 
    wk = sf_complexalloc3 (nh2,ny, nx);   /* k wavefield */
    wx = sf_complexalloc3 (nh, ny, nx);   /* x wavefield */

    ms = sf_intalloc3     (nh, ny, nx);   /* MRS map */
    mr = sf_intalloc3     (nh, ny, nx);   /* MRS map */
    
    skip = sf_boolalloc3 (nrmax,nrmax,nz);/* skip S-R reference slowness combination */

    /* precompute wavenumbers */
    for (ix=0; ix<nx; ix++) {
	x = x01 + ix*dx1;

	iv = 0.5+(x-v0)/dv;
	if      (iv <   0) iv=0;
	else if (iv >= nv) iv=nv-1;
	ii[ix] = iv;
    }

    for (iy=0; iy<ny; iy++) {
	jy = (iy < ny/2)? iy + ny/2: iy - ny/2;
	ky = y0 + jy*dy;
	y = y01 + iy*dy1;
	
	iu = 0.5+(y-u0)/du;
	if      (iu <   0) iu=0;
	else if (iu >= nu) iu=nu-1;
	ij[iy] = iu;

	for (ih=0; ih<nh; ih++) {
	    h = h01 + ih*dh1;
	    
	    iu = 0.5+(y-h-u0)/du;
	    if      (iu <   0) iu=0;
	    else if (iu >= nu) iu=nu-1;
	    is[iy][ih] = iu;
	    
	    iu = 0.5+(y+h-u0)/du;
	    if      (iu <   0) iu=0;
	    else if (iu >= nu) iu=nu-1;
	    ir[iy][ih] = iu;
	}

	for (ih=0; ih<nh2; ih++) {
	    jh = (ih < nh2/2)? ih + nh2/2: ih - nh2/2;
	    kh = h0 + jh*dh;

	    k = 0.5*(ky-kh);
	    ks[iy][ih] = k*k;

	    k = 0.5*(ky+kh);
	    kr[iy][ih] = k*k;
	}
    }    

    /* precompute taper array */
    taper3_init(ntx,nty,nth);

    mms = fslice_init(nh,nx*ny,nz,sizeof(int));
    mmr = fslice_init(nh,nx*ny,nz,sizeof(int));    
}

void comaz_close(void)
/*< free allocated storage >*/
{
    free(**pk); 
    free( *pk); free(pk);
    free(**wk); 
    free( *wk); free(wk);
    free(**wx); 
    free( *wx); free(wx);

    free( *sz); free(sz);
    free( *sm); free(sm);
    free(  nr);

    free(**qq); 
    free( *qq); free( qq);

    free( *ss); free( ss);
    free( *ks); free( ks);
    free( *kr); free( kr);
    free( *is); free( is);
    free( *ir); free( ir);

    fslice_close(mms);
    fslice_close(mmr);
    taper3_close();
    fft3_close();
}

void comaz(bool verb                   /* verbosity flag */, 
	   bool inv                    /* migration/modeling flag */, 
	   float eps                   /* stability factor */,  
	   int nw, float dw, float w0  /* frequency (radian) */,
	   sf_file data                /* data       [nw][nx][ny][nh] */,
	   slice imag                  /* image file [nz][nv][nu][nh] */,
	   slice slow                  /* slowness   [nz][nv][nu]     */,
	   float dt                    /* time error */)
/*< Apply migration/modeling >*/
{
    int iz,iw,ix,iy,ih,iu,iv, j,k;
    int       jx;
    float sy, kx;
    float complex cshift, cref, w, w2, cs, cr, kh, kss, krr;

    if (!inv) { /* prepare image for migration */
	LOOPvuh( qq[iv][iu][ih] = 0.0; );
	for (iz=0; iz<nz; iz++) {
	    slice_put(imag,iz,qq[0][0]);
	}
    }

    /* compute reference slowness */
    for (iz=0; iz<nz; iz++) {
	slice_get(slow,iz,ss[0]);
	nr[iz] = slowref(nrmax,dt/dz,nu*nv,ss[0],sz[iz],sm[iz]);
	if (verb) sf_warning("nr[%d]=%d",iz,nr[iz]);
	/* create MRS map */
	LOOPxyh( ms[ix][iy][ih] = 0; mr[ix][iy][ih] = 0.;);
	for (j=0; j<nr[iz]-1; j++) {
	    sy = 0.5*(sz[iz][j]+sz[iz][j+1]);
	    LOOPxyh( if(ss[ii[ix]][ is[iy][ih] ] > sy) ms[ix][iy][ih]++;
		     if(ss[ii[ix]][ ir[iy][ih] ] > sy) mr[ix][iy][ih]++; );
	}
	fslice_put(mms,iz,ms[0][0]);
	fslice_put(mmr,iz,mr[0][0]);
	for (j=0; j<nr[iz]; j++) {
	    for (k=0; k<nr[iz]; k++) {
		skip[iz][j][k] = true;
	    }
	}
	LOOPxyh( skip[iz][ms[ix][iy][ih]][mr[ix][iy][ih]] = false; );
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
	    slice_get(slow,nz-1,ss[0]);

	    /* imaging condition */
	    slice_get(imag,nz-1,qq[0][0]);	    
	    LOOPxyh( wx[ix][iy][ih] = qq[ ii[ix] ][ ij[iy] ][ih];  );

	    /* loop over migrated depths z */
	    for (iz=nz-2; iz>=0; iz--) {
		/* w-x @ bottom */
		LOOPxyh2( pk[ix][iy][ih] = 0.;);
		taper3(true,true,false,nx,ny,nh,wx);
		LOOPxyh( sy = 0.5*(ss[ ii[ix] ][ is[iy][ih] ] + 
				   ss[ ii[ix] ][ ir[iy][ih] ]);
			 cshift = cexpf(-w*sy*dz);
			 pk[ix][iy][ih] = 
			 wx[ix][iy][ih] * cshift; );

		/* FFT */
		fft3(false,pk);		
		
		slice_get(slow,iz,ss[0]);
		fslice_get(mms,iz,ms[0][0]);
		fslice_get(mmr,iz,mr[0][0]);

		LOOPxyh( wx[ix][iy][ih] = 0; );
		for (j=0; j<nr[iz]; j++) {
		    for (k=0; k<nr[iz]; k++) {
			if (skip[iz][j][k]) continue; /* skip S-R reference combinations */

			/* w-k phase shift */
			cref =        csqrtf(w2*sm[iz][j])
			    +         csqrtf(w2*sm[iz][k]);
			LOOPxyh2(jx = (ix < nx/2)? ix + nx/2: ix - nx/2;
				 kx = x0 + jx*dx; 
				 cs = csqrtf(w2*sm[iz][j]+ks[iy][ih]);
				 cr = csqrtf(w2*sm[iz][k]+kr[iy][ih]);
				 /* comaz approximation */
				 kh = kx*(cr-cs)/(cr+cs); 
				 kss = 0.5*(kx-kh);
				 krr = 0.5*(kx+kh);
				 kss = kss*kss + ks[iy][ih];
				 krr = krr*krr + kr[iy][ih];
				 cs = csqrtf(w2*sm[iz][j] + kss);
				 cr = csqrtf(w2*sm[iz][k] + krr);
				 cshift = cexpf((cref-cs-cr)*dz); 
				 wk[ix][iy][ih] = 
				 pk[ix][iy][ih]*cshift; ); 
			
			/* IFT */
			fft3(true,wk);

			/* accumulate wavefield */
			LOOPxyh( if (ms[ix][iy][ih]==j && 
				     mr[ix][iy][ih]==k) 
				 wx[ix][iy][ih] += wk[ix][iy][ih]; );
		    }
		}
		
		slice_get(imag,iz,qq[0][0]);

		/* w-x at top */
		LOOPxyh( sy = 0.5*(ss[ ii[ix] ][ is[iy][ih] ] + 
				   ss[ ii[ix] ][ ir[iy][ih] ]);
			 cshift = cexpf(-w*sy*dz);
			 wx   [ix]    [iy] [ih] = 
			 qq[ii[ix]][ij[iy]][ih] + 
			 wx   [ix]    [iy] [ih]*cshift; );

	    } /* iz */
	    

	    taper3(true,true,false,nx,ny,nh,wx);
	    sf_complexwrite(wx[0][0],nx*ny*nh,data);

	} else { /* MIGRATION */
	    slice_get(slow,0,ss[0]);

	    sf_complexread(wx[0][0],nx*ny*nh,data);
	    taper3(true,true,false,nx,ny,nh,wx);

	    /* loop over migrated depths z */
	    for (iz=0; iz< nz-1; iz++) {

		/* imaging condition */
		slice_get(imag,iz,qq[0][0]);
		LOOPxyh(        qq[ii[ix]][ij[iy]][ih] += 
			 crealf(wx   [ix]    [iy] [ih] ); );
		slice_put(imag,iz,qq[0][0]);

		/* w-x @ top */
		LOOPxyh2( pk[ix][iy][ih] = 0.;);
		LOOPxyh( sy = 0.5*(ss[ ii[ix] ][ is[iy][ih] ] + 
				   ss[ ii[ix] ][ ir[iy][ih] ]);
			 cshift = conjf(cexpf(-w*sy*dz));
			 pk[ix][iy][ih] = 
			 wx[ix][iy][ih] * cshift; );

		/* FFT */
		fft3(false,pk);

		slice_get(slow,iz+1,ss[0]);
		fslice_get(mms,iz,ms[0][0]);
		fslice_get(mmr,iz,mr[0][0]);
		
		LOOPxyh( wx[ix][iy][ih] = 0; );
		for (j=0; j<nr[iz]; j++) {
		    for (k=0; k<nr[iz]; k++) {
			if (skip[iz][j][k]) continue;
		
			/* w-k phase shift */
			cref = csqrtf(w2*sm[iz][j])
			    +  csqrtf(w2*sm[iz][k]);
			LOOPxyh2(jx = (ix < nx/2)? ix + nx/2: ix - nx/2;
				 kx = x0 + jx*dx; 
				 cs = csqrtf(w2*sm[iz][j] + ks[iy][ih]);
				 cr = csqrtf(w2*sm[iz][k] + kr[iy][ih]);
				 /* comaz approximation */
				 kh = kx*(cr-cs)/(cr+cs); 
				 kss = 0.5*(kx-kh);
				 krr = 0.5*(kx+kh);
				 kss = kss*kss + ks[iy][ih];
				 krr = krr*krr + kr[iy][ih];
				 cs = csqrtf(w2*sm[iz][j] + kss);
				 cr = csqrtf(w2*sm[iz][k] + krr);
				 cshift = conjf(cexpf((cref-cs-cr)*dz)); 
				 wk[ix][iy][ih] = 
				 pk[ix][iy][ih] * cshift; ); 
			
			/* IFT */
			fft3(true,wk);

			/* accumulate wavefield */
			LOOPxyh( if (ms[ix][iy][ih]==j && 
				     mr[ix][iy][ih]==k) 
				 wx[ix][iy][ih] += wk[ix][iy][ih]; );
		    }
		} /* j loop */


		/* w-x @ bottom */
		LOOPxyh( sy = 0.5*(ss[ ii[ix] ][ is[iy][ih] ] + 
				   ss[ ii[ix] ][ ir[iy][ih] ]);
			 cshift = conjf(cexpf(-w*sy*dz));
			 wx[ix][iy][ih] *= cshift; );
		taper3(true,true,false,nx,ny,nh,wx);
	    } /* iz */
	    	    
	    /* imaging condition @ bottom */
	    slice_get(imag,nz-1,qq[0][0]);
	    LOOPxyh(        qq[ii[ix]][ij[iy]][ih] += 
		     crealf(wx   [ix]    [iy] [ih]); );
	    slice_put(imag,nz-1,qq[0][0]);
	} /* else */
    } /* iw */
}

