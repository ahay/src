/* 1-D complex lowrank FFT wave extrapolation */
/*
  Copyright (C) 2008 University of Texas at Austin
  
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
#include <rsf.h>

#include "fft1.h"

int main(int argc, char* argv[]) 
{
    int nx, nx2, nk, nt, m, ix, ik, it, im, n2;
    float dt, f;
    float *curr, **wave;
    sf_complex **lft, **rht, *cwave, *cwavem, **mat;
    sf_file inp, out, prop, left, right;

    sf_init(argc,argv);
    inp = sf_input("in");
    out = sf_output("out");

    if (SF_FLOAT != sf_gettype(inp)) sf_error("Need float input");
    
    if (!sf_histint(inp,"n1",&nx)) sf_error("No n1= in input");

    if (!sf_getint("nt",&nt)) sf_error("No nt= in input");
    if (!sf_getfloat("dt",&dt)) sf_error("No dt= in input");

    sf_putint(out,"n2",nt);
    sf_putfloat(out,"d2",dt);
    sf_putfloat(out,"o2",0.);
    sf_putstring(out,"label2","Time");
    sf_putstring(out,"unit2","s");

    nk = fft1_init(nx,&nx2);

    curr = sf_floatalloc(nx2);
    cwave = sf_complexalloc(nk);

    if (NULL != sf_getstring("right")) {
	left = sf_input("left");   /* Left matrix */
	right = sf_input("right"); /* Right matrix */
	
	if (SF_COMPLEX != sf_gettype(left) ||
	    SF_COMPLEX != sf_gettype(right)) sf_error("Need complex left and right");
	
	if (!sf_histint(left,"n1",&n2) || n2 != nx) sf_error("Need n1=%d in left",nk);
	if (!sf_histint(left,"n2",&m)) sf_error("Need n2=%d in left",m);
	
	if (!sf_histint(right,"n1",&n2) || n2 != m) sf_error("Need n1=%d in right",m);
	if (!sf_histint(right,"n2",&n2) || n2 != nk) sf_error("Need n2=%d in right",nx);
	
	lft = sf_complexalloc2(nx,m);
	rht = sf_complexalloc2(m,nk);
	
	sf_complexread(lft[0],nx*m,left);
	sf_complexread(rht[0],nk*m,right);
	
	sf_fileclose(left);
	sf_fileclose(right);

	cwavem = sf_complexalloc(nk);
	wave = sf_floatalloc2(nx,m);

	ifft1_allocate(cwavem);

	mat = NULL;
    } else { /* Use propagator matrix (for testing) */
	 prop = sf_input("prop");
	 if (SF_COMPLEX != sf_gettype(prop)) sf_error("Need complex prop");

	lft = NULL;
	rht = NULL;

	mat = sf_complexalloc2(nk,nx);
	sf_complexread(mat[0],nx*nk,prop);
	sf_fileclose(prop);

	cwavem = NULL;
	wave = NULL;
    }
 
    /* read the initial data */
    sf_floatread (curr,nx,inp);

    for (ix = nx; ix < nx2; ix++) {
	curr[ix] = 0.;
    }

    /* propagation in time */
    for (it=0; it < nt; it++) {
	/* FFT: curr -> cwave */
	fft1(curr,cwave);

	if (NULL == mat) {
	    for (im = 0; im < m; im++) {
		for (ik = 0; ik < nk; ik++) {
#ifdef SF_HAS_COMPLEX_H
		    cwavem[ik] = cwave[ik]*rht[ik][im];
#else
		    cwavem[ik] = sf_cmul(cwave[ik],rht[ik][im]);
#endif
		}
		/* Inverse FFT: cwavem -> wave[im] */
		ifft1(wave[im],cwavem);
	    }
	}

	for (ix = 0; ix < nx; ix++) {
	    f = 0.0f;
	    if (NULL == mat) {
		for (im = 0; im < m; im++) {
		    f += crealf(lft[im][ix])*wave[im][ix];
		}
	    } else {
		for (ik = 0; ik < nk; ik++) {
#ifdef SF_HAS_COMPLEX_H
		    f += crealf(mat[ix][ik] * cwave[ik]);
#else
		    f += crealf(sf_cmul(mat[ix][ik],cwave[ik]));
#endif
		}
	    }
	    curr[ix] = f;
	} 

	sf_floatwrite(curr,nx,out);
    }
    
    exit(0);
}
