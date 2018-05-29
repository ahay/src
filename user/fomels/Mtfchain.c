/* Find a chain of Fourier weighting and scaling */
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

#include "fchain.h"
#include "twosmooth.h"

int main(int argc, char* argv[])
{
    int i, n, nw, n2, rect, frect, iter, niter, liter;
    float *w, *dw, *x, *y, *r, *p;
    sf_file wht, fwht, src, tgt, mch;

    sf_init(argc,argv);
    src = sf_input("in");
    wht = sf_output("out");

    tgt = sf_input("target");
    fwht = sf_output("fweight");
    mch = sf_output("match");

    if (!sf_histint(src,"n1",&n)) sf_error("No n1= in input");

    nw = kiss_fft_next_fast_size((n+1)/2)+1;

    n2 = n+nw;
    sf_putint(fwht,"n1",nw);

    w = sf_floatalloc(n2);
    dw = sf_floatalloc(n2);

    x = sf_floatalloc(n);
    y = sf_floatalloc(n);
    r = sf_floatalloc(n);

    if (!sf_getint("rect",&rect)) rect=1;
    /* smoothing in time */
    if (!sf_getint("frect",&frect)) frect=1;
    /* smoothing in frequency */

    twosmooth_init(n,nw,rect,frect);

    sf_floatread(x,n,src);
    sf_floatread(y,n,tgt);

    fchain_init(n,nw,w,w+n,x);

    sf_conjgrad_init(n2, n2, n, n, 1., 1.e-6, true, false);

    p = sf_floatalloc(n2);

    /* initialize */
    for (i=0; i < n2; i++) {
	w[i] = 1.0f;
    }

    if (!sf_getint("niter",&niter)) niter=0;
    /* number of iterations */
    if (!sf_getint("liter",&liter)) liter=50;
    /* number of linear iterations */

    for (iter=0; iter < niter; iter++) {
	fchain_apply(r);
	
	for (i=0; i < n; i++) {
	    r[i] = y[i] - r[i];
	}
	
	sf_conjgrad(NULL, fchain_lop,twosmooth_lop,p,dw,r,liter);
	
	for (i=0; i < n2; i++) {
	    w[i] += dw[i];
	}
    }

    sf_floatwrite(w,n,wht);
    sf_floatwrite(w+n,nw,fwht);
    
    fchain_apply(r);
    sf_floatwrite(r,n,mch);

    exit(0);
}
