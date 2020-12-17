/* Find a symmetric chain of Fourier weighting and scaling */
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

#include <rsf.h>

#include "sfchain.h"
#include "twosmooth.h"

int main(int argc, char* argv[])
{
    int i, n, nw, n2, rect, frect, iter, niter, liter;
    float *w, *dw, *x, *y, *r, *p; /* *num, *den; */
    sf_file wht, fwht, src, tgt, mch;

    sf_init(argc,argv);
    src = sf_input("in");
    wht = sf_output("out");

    tgt = sf_input("target");
    fwht = sf_output("fweight");
    mch = sf_output("match");

    if (!sf_histint(src,"n1",&n)) sf_error("No n1= in input");

    nw = kiss_fft_next_fast_size((n+1)/2)+1;

    n2 = 3*n+nw;
    sf_putint(fwht,"n1",nw);

    w = sf_floatalloc(n2);
    dw = sf_floatalloc(n2);

    x = sf_floatalloc(n);
    y = sf_floatalloc(n);
    r = sf_floatalloc(3*n);

    /*
    num = sf_floatalloc(nw);
    den = sf_floatalloc(nw);
    */

    if (!sf_getint("rect",&rect)) rect=1;
    /* smoothing in time */
    if (!sf_getint("frect",&frect)) frect=1;
    /* smoothing in frequency */

    sf_floatread(x,n,src);
    sf_floatread(y,n,tgt);

    sfchain_init(n,nw,w+2*n,w+3*n,w,w+n,x);

    p = sf_floatalloc(n2);

    /* initialize */
    for (i=0; i < 2*n; i++) {
	w[i] = 0.0f;
    }
    for (i=2*n; i < n2; i++) {
	w[i] = 1.0f;
    }

    if (!sf_getint("niter",&niter)) niter=0;
    /* number of iterations */
    if (!sf_getint("liter",&liter)) liter=50;
    /* number of linear iterations */

    /* smart initialization */
    
    /*    sf_divn_init(1, n, &n, &rect, liter, true);
    sf_divne (y, x, w+2*n, 0.0f);
    for (i=2*n; i < 3*n; i++) {
	w[i] = sqrtf(fabsf(w[i]));
    }
    sf_divn_close(); */

    /*
    sfchain_div(y,num,den);
    sf_divn_init(1, nw, &nw, &frect, liter, true);
    sf_divne (num, den, w+3*n, 0.0f);
    sf_divn_close(); */
    
    sf_conjgrad_init(n2, n2, 3*n, 3*n, 1., 1.e-6, true, false);
    twosmooth_init(n,nw,rect,frect,2*n);
    
    for (iter=0; iter < niter; iter++) {
	sfchain_res(y,r);
	
	sf_conjgrad(NULL, sfchain_lop,twosmooth_lop,p,dw,r,liter);
	
	for (i=0; i < n2; i++) {
	    w[i] += dw[i];
	}

	/* variable projection 
	sfchain_res(y,r);
	sf_solver(sfchainx_lop,sf_cgstep,2*n,3*n,dw,r,liter,"verb",true,"end");
	sf_cgstep_close();
	for (i=0; i < 2*n; i++) {
	    w[i] += dw[i];
	    } */
    }

    sf_floatwrite(w+2*n,n,wht);
    sf_floatwrite(w+3*n,nw,fwht);
    
    sfchain_apply(y);
    sf_floatwrite(y,n,mch);

    exit(0);
}
