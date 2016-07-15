/* Multi-dimensional PEF (prediction error filter) estimation on a helix. */
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

#include "bound.h"
#include "misinput.h"
#include "pef.h"
#include "compress.h"
#include "printfilter.h"
#include "createhelix.h" 

int main(int argc, char* argv[])
{
    int n[SF_MAX_DIM], n0[SF_MAX_DIM];
    int a[SF_MAX_DIM], center[SF_MAX_DIM], gap[SF_MAX_DIM];
    int dim, n123, i, niter, na, *kk;
    float *dd, tol;
    sf_filter aa;   
    char varname[6], *lagfile;
    sf_file in, filt, lag, mask;

    sf_init (argc, argv);
    in = sf_input("in");
    filt = sf_output("out");

    if (NULL == (lagfile = sf_getstring("lag"))) sf_error("Need lag=");
    /* output file for filter lags */

    lag = sf_output(lagfile);
    sf_settype(lag,SF_INT);

    sf_putstring(filt,"lag",lagfile);

    dim = sf_filedims(in,n);
    sf_putints (lag,"n",n,dim);

    if (!sf_getints("a",a,dim)) sf_error("Need a=");

    if (!sf_getints("center",center,dim)) {
	for (i=0; i < dim; i++) {
	    center[i] = (i+1 < dim && a[i+1] > 1)? a[i]/2: 0;
	}
    }

    if (!sf_getint("na",&na)) na=0;
    /* filter size */

    if (!sf_getfloat("tol",&tol)) tol=1.e-6;
    /* tolerance for filter compression */

    if (0 == na) {
	if (!sf_getints("gap",gap,dim)) {
	    for (i=0; i < dim; i++) {
		gap[i] = 0;
	    }
	}

	aa = createhelix(dim, n, center, gap, a); /* allocate PEF */

	for (i=0; i < dim; i++) {	    
	    n0[i] = n[i];
	}
    } else {
	aa =  sf_allocatehelix (na);
	if (!sf_getints ("lags", aa->lag, na)) sf_error("Need lags=");
	if (!sf_getints ("n", n0, dim)) {
	    for (i=0; i < dim; i++) {	    
		n0[i] = n[i];
	    }
	}
    }

    n123 = 1;
    for (i=0; i < dim; i++) {
	n123 *= n[i];
    }

    dd = sf_floatalloc(n123);
    kk = sf_intalloc(n123);

    if (NULL != sf_getstring("maskin")) {
	/* optional input mask file */
	mask = sf_input("maskin");

	switch (sf_gettype(mask)) {
	    case SF_INT:
		sf_intread (kk,n123,mask);
		break;
	    case SF_FLOAT:
		sf_floatread (dd,n123,mask);
		for (i=0; i < n123; i++) {
		    kk[i] = (dd[i] != 0.);
		}
		break;
	    default:
		sf_error ("Wrong data type in maskin");
		break;
	}

	sf_fileclose (mask);
    } else {
	for (i=0; i < n123; i++) {
	    kk[i] = 1;
	}
    }

    sf_floatread (dd,n123,in);

    bound (dim, n0, n, a, aa); 
 	
    find_mask(n123, kk, aa);   /* account for missing data */
    
    if (NULL != sf_getstring("maskout")) {
	/* optional output mask file */
	mask = sf_output("maskout");

	for (i=0; i < n123; i++) {
	    kk[i] = aa->mis[i]? 0: 1;
	}
	
	sf_settype(mask,SF_INT);
	sf_intwrite (kk,n123,mask);
    }

    if(!sf_getint("niter",&niter)) niter=2*(aa->nh);
    /* number of iterations */

    find_pef (n123, dd, aa, niter);         /* estimate aa */
    aa = compress( aa, tol);              /* eliminate zeroes */
    print(dim, n, center, a, aa);           /* print filter */

    sf_putint(filt,"n1",aa->nh);
    sf_putint(lag,"n1",aa->nh);

    for (i=1; i < dim; i++) {
	sprintf(varname,"n%d",i+1);
	sf_putint(filt,varname,1);
	sf_putint(lag,varname,1);
    }

    sf_intwrite(aa->lag,aa->nh,lag);
    sf_floatwrite(aa->flt,aa->nh,filt);


    exit (0);
}

/* 	$Id$	 */
