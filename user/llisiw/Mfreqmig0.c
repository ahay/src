/* Image interpolation coefficients estimation (single-arrival) */
/*
  Copyright (C) 2014 University of Texas at Austin
  
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

#include "freqmig0.h"

int main(int argc, char* argv[])
{
    bool verb, adj, hessian;
    int n1, n2, nf, i;
    float df, of;
    sf_complex **img, *coe0, *coe, *dx, *dr;
    int iter, niter, cgiter, count;
    float rhsnorm, rhsnorm0, rhsnorm1, rate, gamma;
    char *what;
    sf_file in, out, image, coef, grad;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (NULL == (what = sf_getstring("what"))) what="tomo";
    /* what to compute (default tomography) */

    switch (what[0]) {
	case 'l': /* linear operator */

	    if (!sf_getbool("adj",&adj)) adj=false;
	    /* adjoint flag (for what=linear) */

	    /* read image */
	    if (NULL == sf_getstring("image"))
		sf_error("Need image image=");
	    image = sf_input("image");

	    if (!sf_histint(image,"n1",&n1)) sf_error("No n1 in image.");
	    if (!sf_histint(image,"n2",&n2)) sf_error("No n2 in image.");

	    if (!sf_histint(image,"n3",&nf)) sf_error("No n3 in image.");
	    if (!sf_histfloat(image,"o3",&of)) sf_error("No o3 in image.");
	    if (!sf_histfloat(image,"d3",&df)) sf_error("No d3 in image.");

	    img = NULL;

	    /* read coefficients */
	    if (NULL == sf_getstring("coef"))
		sf_error("Need coefficients coef=");
	    coef = sf_input("coef");

	    coe0 = sf_complexalloc(n1*n2*2);
	    sf_complexread(coe0,n1*n2*2,coef);
	    sf_fileclose(coef);
	    
	    /* allocate temporary memory */
	    dx  = sf_complexalloc(n1*n2*2);
	    dr  = sf_complexalloc(n1*n2*nf);

	    /* read input and write output header */
	    if (adj) {
		sf_complexread(dr,n1*n2*nf,in);

		sf_putint(out,"n3",2);
	    } else {
		sf_complexread(dx,n1*n2*2,in);

		sf_putint(out,"n3",nf);
	    }

	    /* initialize operator */
	    freqmig0_init(n1,n2, of,df,nf, img);

	    /* set operator */
	    freqmig0_set(coe0);

	    /* linear operator */
	    if (adj) {
		freqmig0_oper(true,false,n1*n2*2,n1*n2*nf,dx,dr);
	    } else {
		freqmig0_oper(false,false,n1*n2*2,n1*n2*nf,dx,dr);
	    }
	    	    
	    /* write output */
	    if (adj) {
		sf_complexwrite(dx,n1*n2*2,out);
	    } else {
		sf_complexwrite(dr,n1*n2*nf,out);
	    }

	    break;

	case 't': /* tomography */

	    /* read model dimension */
	    if (!sf_histint(in,"n1",&n1)) sf_error("No n1 in input.");
	    if (!sf_histint(in,"n2",&n2)) sf_error("No n2 in input.");
	    
	    /* read initial guess */
	    coe0 = sf_complexalloc(n1*n2*2);
	    sf_complexread(coe0,n1*n2*2,in);
	    
	    /* read image */
	    if (NULL == sf_getstring("image"))
		sf_error("Need image image=");
	    image = sf_input("image");

	    if (!sf_histint(image,"n3",&nf)) sf_error("No n3 in image.");
	    if (!sf_histfloat(image,"o3",&of)) sf_error("No o3 in image.");
	    if (!sf_histfloat(image,"d3",&df)) sf_error("No d3 in image.");

	    img = sf_complexalloc2(n1*n2,nf);
	    sf_complexread(img[0],n1*n2*nf,image);
	    sf_fileclose(image);

	    /* allocate temporary memory */
	    dx  = sf_complexalloc(n1*n2*2);
	    dr  = sf_complexalloc(n1*n2*nf);
	    coe = sf_complexalloc(n1*n2*2);

	    if (!sf_getint("niter",&niter)) niter=5;
	    /* number of inversion iterations */
	    
	    if (!sf_getint("cgiter",&cgiter)) cgiter=10;
	    /* number of conjugate-gradient iterations */

	    if (!sf_getbool("hessian",&hessian)) hessian=false;
	    /* if y, exact Hessian; n, Gauss-Newton */

	    if (!sf_getbool("verb",&verb)) verb=false;
	    /* verbosity flag */

	    /* output gradient at each iteration */
	    if (NULL != sf_getstring("grad")) {
		grad = sf_output("grad");
		sf_putint(grad,"n4",niter);
	    } else {
		grad = NULL;
	    }

	    /* initialize */
	    freqmig0_init(n1,n2, of,df,nf, img);

	    rhsnorm0 = freqmig0_cost(coe0,dr);
	    rhsnorm = rhsnorm0;
	    rhsnorm1 = rhsnorm;
	    rate = rhsnorm1/rhsnorm0;

	    sf_warning("L2 misfit after iteration 0 of %d: %g",niter,rate);

	    /* iterations over inversion */
	    for (iter=0; iter < niter; iter++) {

		/* set operator */
		freqmig0_set(coe0);

		/* solve update */
		sf_csolver(freqmig0_oper,sf_ccgstep,n1*n2*2,n1*n2*nf,dx,dr,cgiter,"verb",verb,"end");
		sf_ccgstep_close();

		/* output gradient */
		if (grad != NULL) sf_complexwrite(dx,n1*n2*2,grad);

		/* line search */
		gamma = 1.;
		for (count=0; count < 5; count++) {
	    
		    /* update */
		    for (i=0; i < n1*n2; i++) {
			coe[i]       = coe0[i]+gamma*sf_cmplx(crealf(dx[i]),0.);
			coe[n1*n2+i] = coe0[n1*n2+i]+gamma*sf_cmplx(crealf(dx[n1*n2+i]),0.);
		    }

		    /* check cost */
		    rhsnorm = freqmig0_cost(coe,dr);
		    rate = rhsnorm/rhsnorm1;
	    
		    if (rate < 1.) {
			for (i=0; i < n1*n2; i++) {
			    coe0[i] = coe[i];
			    coe0[n1*n2+i] = coe[n1*n2+i];
			}
		
			rhsnorm1 = rhsnorm;
			rate = rhsnorm1/rhsnorm0;
			break;
		    }

		    gamma *= 0.5;
		}

		if (count == 5) {
		    sf_warning("Line-search failure at iteration %d of %d.",iter+1,niter);
		    break;
		}

		sf_warning("L2 misfit after iteration %d of %d: %g (line-search %d)",iter+1,niter,rate,count);
	    }

	    /* write output */
	    sf_complexwrite(coe0,n1*n2*2,out);

	    /* clean-up */
	    freqmig0_close();
	    free(dx); free(dr);

	    break;
    }

    exit(0);
}
