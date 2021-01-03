/* 1-D complex lowrank FFT wave extrapolation using complex to complex fft using initial condition*/
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

#include "cfft1.h"

int main(int argc, char* argv[]) 
{
    bool sub;
    int nx, nx2, nk, nt, m, ix, ik, it, im, n2;
    sf_complex *curr, **wave, c;
    float *rcurr, dt;
    sf_complex **lft, **rht, *cwave, *cwavem, **mat;
    sf_file Fw, Fo, right, left, prop;

    sf_init(argc,argv);
    Fw = sf_input("in");
    Fo = sf_output("out");

    if (SF_COMPLEX != sf_gettype(Fw)) sf_error("Need complex input");

    sf_settype(Fo,SF_FLOAT);

    /* Read/Write axes */

    if (!sf_getbool("sub",&sub)) sub=false;
    /* if -1 is included in the matrix */

    if (!sf_histint(Fw,"n1",&nx)) sf_error("No n1= in input");

    if (!sf_getint("nt",&nt)) sf_error("No nt= in input");
    if (!sf_getfloat("dt",&dt)) sf_error("No dt= in input");

    sf_putint(Fo,"n2",nt);
    sf_putfloat(Fo,"d2",dt);
    sf_putfloat(Fo,"o2",0.);
    sf_putstring(Fo,"label2","Time");
    sf_putstring(Fo,"unit2","s");

    nk = cfft1_init(nx,&nx2);
    sf_warning("nk=%d\n", nk);

    curr = sf_complexalloc(nx2);
    rcurr= sf_floatalloc(nx2);
    cwave = sf_complexalloc(nk);
 
    if (NULL != sf_getstring("right")) {
	left = sf_input("left");   /* Left matrix */
	right = sf_input("right"); /* Right matrix */
	
	if (SF_COMPLEX != sf_gettype(left) ||
	    SF_COMPLEX != sf_gettype(right)) sf_error("Need complex left and right");
	
	if (!sf_histint(left,"n1",&n2) || n2 != nx) sf_error("Need n1=%d in left",nx);
	if (!sf_histint(left,"n2",&m)) sf_error("Need n2= in left");
	
	if (!sf_histint(right,"n1",&n2) || n2 != m) sf_error("Need n1=%d in right",m);
	if (!sf_histint(right,"n2",&n2) || n2 != nk) sf_error("Need n2=%d in right, now n2=%d",nk,n2);
	
	lft = sf_complexalloc2(nx,m);
	rht = sf_complexalloc2(m,nk);
	
	sf_complexread(lft[0],nx*m,left);
	sf_complexread(rht[0],nk*m,right);
	
	sf_fileclose(left);
	sf_fileclose(right);
	
	cwavem = sf_complexalloc(nk);
	wave = sf_complexalloc2(nx2,m);

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
	
    sf_complexread (curr,nx,Fw);

    for (ix = nx; ix < nx2; ix++) {
	curr[ix] = sf_cmplx(0.,0.);
    }

    for (ix = 0; ix < nx2; ix++) {
	rcurr[ix]=0.;
    }

    /* propagation in time */
    for (it=0; it < nt; it++) {
	sf_warning("it=%d;",it);

        /* FFT: curr -> cwave */
	cfft1(curr,cwave);

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
		icfft1(wave[im],cwavem);
	    }
	}

	for (ix = 0; ix < nx; ix++) {

	    c = sub? curr[ix]: sf_cmplx(0.,0.);

	    if (NULL == mat) {
		for (im = 0; im < m; im++) {
#ifdef SF_HAS_COMPLEX_H
		    c += lft[im][ix]*wave[im][ix];
#else
		    c = sf_cadd(c,sf_cmul(lft[im][ix], wave[im][ix]));
#endif
		}
	    } else {
		for (ik = 0; ik < nk; ik++) {
#ifdef SF_HAS_COMPLEX_H
		    c += mat[ix][ik] * cwave[ik];
#else
		    c = sf_cadd(c,sf_cmul(mat[ix][ik],cwave[ik]));
#endif
		}
	    }

	    curr[ix] = c;
	    rcurr[ix] = crealf(c);
	} 

	sf_floatwrite(rcurr,nx,Fo);
    }
    
    exit(0);
}
