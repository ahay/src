/* Find a symmetric chain of 2D-Fourier weighting and scaling */
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
#include "chain2dfft.h"
#include "twosmooth2.h"
#include "fft2.h"

int main(int argc, char* argv[])
{
    int i, n, nk, n2, iter, niter, liter, nt, nx, nt1, nt2, nx2;
    int rect1, rect2, frect1, frect2; 
    float *w, *dw, *x, *y, *r, *p;
    sf_file wht, fwht, src, tgt, mch;
    /* For fft2 */
    bool isCmplx = false;
    int pad = 1;

    sf_init(argc,argv);
    src = sf_input("in");
    wht = sf_output("out"); /* space weight */

    tgt = sf_input("target");
    fwht = sf_output("fweight"); /* frequency weight */
    mch = sf_output("match"); /* matched */

    if (!sf_histint(src,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histint(src,"n2",&nx)) sf_error("No n2= in input");
    n = nt*nx; 

    nk = fft2_init(isCmplx, pad, nt, nx, &nt1, &nx2);
    nt2 = nk/nx2;
 
    n2 = 3*n+nk; 
    sf_putint(fwht,"n1",nt2);
    sf_putint(fwht,"n2",nx2);

    w = sf_floatalloc(n2); 
    dw = sf_floatalloc(n2);

    x = sf_floatalloc(n); 
    y = sf_floatalloc(n);
    r = sf_floatalloc(3*n);

    if (!sf_getint("rect1",&rect1)) rect1=1;
    if (!sf_getint("rect2",&rect2)) rect2=1;
    /* smoothing in time */
    if (!sf_getint("frect1",&frect1)) frect1=1;
    if (!sf_getint("frect2",&frect2)) frect2=1;
    /* smoothing in frequency */

    twosmooth2_init(n,nk,nt,nt2,
            rect1,rect2,
            frect1,frect2,
            2*n);

    sf_floatread(x,n,src); /* source */
    sf_floatread(y,n,tgt); /* target */

    sfchain2d_init(nt,nx,nt1,nx2,nk,
		   w+2*n,w+3*n,w,w+n,x);

    sf_conjgrad_init(n2, n2, 3*n, 3*n, 1., 1.e-6, true, false);
 
    p = sf_floatalloc(n2); 

    /* initialize w [time w and freqz w] */
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

    for (iter=0; iter < niter; iter++) {
    sf_warning("Start %d",iter);

    sfchain2d_res(y,r);

    sf_warning("Residual %d",iter);
    
    sf_conjgrad(NULL, sfchain2d_lop,twosmooth2_lop,p,dw,r,liter);


    for (i=0; i < n2; i++) {
        w[i] += dw[i];
    }   
    }

    sf_floatwrite(w+2*n,n,wht); 
    sf_floatwrite(w+3*n,nk,fwht);
    sfchain2d_apply(y);
    sf_floatwrite(y,n,mch);

    exit(0);
}
