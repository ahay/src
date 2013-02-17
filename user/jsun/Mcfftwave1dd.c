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

#include "cfft1.h"

int main(int argc, char* argv[]) 
{
    int nx, nx2, nk, nt, m, ix, ik, it, im, n2;
    float dt, fscale;
    sf_complex *curr, **wave,f;
    float *rcurr;
    sf_complex **lft, **rht, *cwave, *cwavem;
    sf_file inp, out, left, right, refl;
    float *rr,*ww;

    sf_init(argc,argv);
    inp = sf_input("in");
    out = sf_output("out");
    refl = sf_input("refl");

    if (SF_COMPLEX != sf_gettype(inp)) sf_error("Need complex input");
    
    if (!sf_histint(refl,"n1",&nx)) sf_error("No n1= in input");

    if (!sf_histint(inp,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histfloat(inp,"d1",&dt)) sf_error("No d1= in input");
    sf_putint(out,"n1",nx);
    sf_putint(out,"n2",nt);
    sf_putfloat(out,"d2",dt);
    sf_putfloat(out,"o2",0.);
    sf_putstring(out,"label2","Time");
    sf_putstring(out,"unit2","s");

    nk = fft1_init(1,nx,&nx2);
    sf_warning("nk=%d\n", nk);
     
    curr = sf_complexalloc(nx2);
    rcurr= sf_floatalloc(nx2);
    cwave = sf_complexalloc(nk);

    left = sf_input("left");   /* Left matrix */
    right = sf_input("right"); /* Right matrix */

    if (SF_COMPLEX != sf_gettype(left) ||
	SF_COMPLEX != sf_gettype(right)) sf_error("Need complex left and right");

    if (!sf_histint(left,"n1",&n2) || n2 != nx) sf_error("Need n1=%d in left",nk);
    if (!sf_histint(left,"n2",&m)) sf_error("Need n2=%d in left",m);

    if (!sf_histint(right,"n1",&n2) || n2 != m) sf_error("Need n1=%d in right",m);
    if (!sf_histint(right,"n2",&n2) || n2 != nk) sf_error("Need n2=%d in right, now n2=%d",nx,n2);

    lft = sf_complexalloc2(nx,m);
    rht = sf_complexalloc2(m,nk);

    sf_complexread(lft[0],nx*m,left);
    sf_complexread(rht[0],nk*m,right);
	
    sf_fileclose(left);
    sf_fileclose(right);

    cwavem = sf_complexalloc(nk);
    wave = sf_complexalloc2(nx,m);

    /* read wavelet & reflectivity */
    ww=sf_floatalloc(nt);  sf_floatread(ww,nt,inp);
    rr=sf_floatalloc(nx); sf_floatread(rr,nx,refl);

    ifft1_allocate(cwavem);

    /*  
    sf_floatread (curr,nx,inp);

    for (ix = nx; ix < nx2; ix++) {
	curr[ix] = 0.;
    }
    */


    for (ix = 0; ix < nx2; ix++) {
	curr[ix] = sf_cmplx(0.,0.);
	rcurr[ix]= 0.;
    }

    fscale = 1.0/nx2;

    /* propagation in time */
    for (it=0; it < nt; it++) {
	/* FFT: curr -> cwave */
	fft1(curr,cwave);

	for (im = 0; im < m; im++) {
	    for (ik = 0; ik < nk; ik++) {
#ifdef SF_HAS_COMPLEX_H
		cwavem[ik] = cwave[ik]*rht[ik][im]*fscale;
#else
		cwavem[ik] = sf_crmul(sf_cmul(cwave[ik],rht[ik][im]),fscale);
#endif
		//sf_warning("realcwave=%g, imagcwave=%g", crealf(cwavem[ik]),cimagf(cwavem[ik]));
	    }
	    /* Inverse FFT: cwavem -> wave[im] */
	    ifft1(wave[im],cwavem);
	}

	for (ix = 0; ix < nx; ix++) {
#ifdef SF_HAS_COMPLEX_H
	    f = ww[it] * rr[ix]; // source term
#else
	    f = sf_crmul(ww[it],rr[ix]);
#endif
	    for (im = 0; im < m; im++) {
#ifdef SF_HAS_COMPLEX_H
		f += lft[im][ix]*wave[im][ix];
#else
		f += sf_cmul(lft[im][ix], wave[im][ix]);
#endif
	    }
	    curr[ix] = f;
	    //sf_warning("f= %g", f);
	    rcurr[ix] = crealf(f);
	} 

	sf_floatwrite(rcurr,nx,out);
    }
    
    exit(0);
}
