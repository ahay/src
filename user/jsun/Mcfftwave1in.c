/* 1-D complex lowrank FFT wave extrapolation using complex to complex fft BY INJECTION*/
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
    int nx, nx2, nk, nt, m, ix, ik, it, im, n2;
    sf_complex *curr, **wave, c, *ww;
    float *rcurr;
    sf_complex **lft, **rht, *cwave, *cwavem;
    sf_file Fw, Fr, Fo, right, left;
    sf_axis at,ax;
    float *rr;

    sf_init(argc,argv);
    Fw = sf_input("in");
    Fo = sf_output("out");
    Fr = sf_input("refl");

    if (SF_COMPLEX != sf_gettype(Fw)) sf_error("Need complex input");
    if (SF_FLOAT != sf_gettype(Fr)) sf_error("Need float refl");

    sf_settype(Fo,SF_FLOAT);

    /* Read/Write axes */
    at = sf_iaxa(Fw,1); nt = sf_n(at);
    ax = sf_iaxa(Fr,1); nx = sf_n(ax);

    sf_oaxa(Fo,ax,1);
    sf_oaxa(Fo,at,2);


    nk = cfft1_init(nx,&nx2);
    sf_warning("nk=%d\n", nk);
     
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

    curr = sf_complexalloc(nx2);
    rcurr= sf_floatalloc(nx2);

    cwave = sf_complexalloc(nk);
    cwavem = sf_complexalloc(nk);
    wave = sf_complexalloc2(nx2,m);

    /* read wavelet & reflectivity */
    ww=sf_complexalloc(nt);  sf_complexread(ww,nt,Fw);
    rr=sf_floatalloc(nx); sf_floatread(rr,nx,Fr);

    for (ix = 0; ix < nx2; ix++) {
	curr[ix] = sf_cmplx(0.,0.);
	rcurr[ix]= 0.;
    }


    /* propagation in time */
    for (it=0; it < nt; it++) {
	sf_warning("it=%d;",it);

        /* FFT: curr -> cwave */
	cfft1(curr,cwave);

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

	for (ix = 0; ix < nx; ix++) {
#ifdef SF_HAS_COMPLEX_H
	    c = ww[it] * rr[ix]; // source term
#else
	    c = sf_crmul(ww[it],rr[ix]);
#endif
	    for (im = 0; im < m; im++) {
#ifdef SF_HAS_COMPLEX_H
		c += lft[im][ix]*wave[im][ix];
#else
                c = sf_cadd(c,sf_cmul(lft[im][ix], wave[im][ix]));
#endif
	    }
	    curr[ix] = c;
	    rcurr[ix] = crealf(c);
//	    rcurr[ix] = cimagf(c);
	} 

	sf_floatwrite(rcurr,nx,Fo);
    }
    
    exit(0);
}
