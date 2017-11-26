/* Multi-dimensional missing data interpolation with shaping regularization. */
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

#include "hshape.h"
#include "bound.h"

int main(int argc, char* argv[])
{
    int i, ia, na, nx, ns, dim, niter;
    int n[SF_MAX_DIM], m[SF_MAX_DIM], a[SF_MAX_DIM];
    float a0, eps, *mm, *pp;
    bool verb, *known;
    sf_filter aa;
    char* lagfile;
    sf_file in, out, filt, lag, mask;

    sf_init (argc,argv);
    in = sf_input("in");
    filt = sf_input("filt");
    out = sf_output("out");

    dim = sf_filedims (in,n);

    if (!sf_histint(filt,"n1",&na)) sf_error("No n1= in filt");
    aa = sf_allocatehelix (na);

    if (!sf_histfloat(filt,"a0",&a0)) a0=1.;
    if (!sf_histints(filt,"a",a,dim)) {
	for (i=0; i < dim; i++) {
	    a[i]=1;
	}
    }

    if (NULL != (lagfile = sf_getstring("lag")) /* file with filter lags */
	|| 
	NULL != (lagfile = sf_histstring(filt,"lag"))) {
	lag = sf_input(lagfile);

	sf_intread(aa->lag,na,lag);
    } else {
	lag = NULL;
	for( ia=0; ia < na; ia++) {
	    aa->lag[ia] = ia+1;
	}
    }

    
    if (!sf_getints ("n",m,dim) && (NULL == lag ||
				    !sf_histints (lag,"n",m,dim))) {
	for (i=0; i < dim; i++) {
	    m[i] = n[i];
	}
    }
 
    if (NULL != lag) sf_fileclose(lag);

    bound (dim, false, m, n, a, aa);

    sf_floatread (aa->flt,na,filt);
    sf_fileclose(filt);
    
    for( ia=0; ia < na; ia++) {
	aa->flt[ia] /= a0;
    }

    if (!sf_getint ("ns",&ns)) sf_error("Need ns=");
    /* scaling */
    if (!sf_getint("niter",&niter)) niter=100;
    /* Number of iterations */
    if (!sf_getfloat("eps",&eps)) eps=1.;
    /* regularization parameter */
    if (!sf_getbool("verb",&verb)) verb=true;
    /* verbosity flag */

    nx = 1;
    for( i=0; i < dim; i++) {
	nx *= n[i];
    }
  
    mm = sf_floatalloc (nx);
    pp = sf_floatalloc (nx);
    known = sf_boolalloc (nx);

    sf_mask_init(known);
    hshape_init (nx,ns,aa);
    sf_conjgrad_init(nx, nx, nx, nx, eps, 1.e-8, verb, false);

    if (NULL != sf_getstring("mask")) {
	/* optional input mask file for known data */
	mask = sf_input("mask");
	sf_floatread(mm,nx,mask);
	sf_fileclose(mask);
	
	for (i=0; i < nx; i++) {
	    known[i] = (bool) (mm[i] != 0.);
	}
	sf_floatread(mm,nx,in);
    } else {
	sf_floatread(mm,nx,in);
	   
	for (i=0; i < nx; i++) {
	    known[i] = (bool) (mm[i] != 0.);
	}
    }

    sf_conjgrad(NULL, sf_mask_lop, hshape_lop, pp, mm, mm, niter);
    
    sf_floatwrite (mm,nx,out);


    exit (0);
}


