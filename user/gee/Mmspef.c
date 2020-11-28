/* Multi-scale PEF estimation.
*/
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

#include "mshelix.h"
#include "createmshelix.h"
#include "mspef.h"
#include "misinput.h"
#include "printfilter.h"

int main(int argc, char* argv[])
{
    int n123, is,ns, dim,i, niter, *jump, *kk, nh;
    int n[SF_MAX_DIM], a[SF_MAX_DIM], center[SF_MAX_DIM], gap[SF_MAX_DIM];
    float *dd;
    msfilter msaa;
    char varname[6], *lagfile;
    sf_file in, pef, lag, mask;

    sf_init(argc,argv);
    in = sf_input("in");
    pef = sf_output("out");

    if (NULL == (lagfile = sf_getstring("lag"))) sf_error("Need lag=");
    /* output file for filter lags */

    lag = sf_output(lagfile);
    sf_settype(lag,SF_INT);

    sf_putstring(pef,"lag",lagfile);

    dim = sf_filedims(in,n);
    sf_putints (lag,"n",n,dim);

    if (!sf_getints("a",a,dim)) sf_error("Need a=");
    sf_putints (pef,"a",a,dim);

    if (!sf_getints("center",center,dim)) {
	for (i=0; i < dim; i++) {
	    center[i] = (i+1 < dim && a[i+1] > 1)? a[i]/2: 0;
	}
    }

    if (!sf_getints("gap",gap,dim)) {
	for (i=0; i < dim; i++) {
	    gap[i] = 0;
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
	    kk[i] = 1.;
	}
    }

    sf_floatread (dd,n123,in);

    if (!sf_getint("ns",&ns)) sf_error("Need ns=");
    /* number of scales */
    jump = sf_intalloc(ns);

    if (!sf_getints("jump",jump,ns)) sf_error("Need jump=");

    msaa = createmshelix(dim, n, center, gap, ns, jump, a); /* allocate PEF */

    nh = msaa->nh;
    if (!sf_getint("niter",&niter)) niter=nh*2;

    for (is=0; is < ns; is++) {
	onescale (is, msaa);
	find_mask(n123, kk, msaa->one); /* missing data */
    }

    if (NULL != sf_getstring("maskout")) {
	/* optional output mask file */
	mask = sf_output("maskout");	
	sf_settype(mask,SF_INT);
	
	sprintf(varname,"n%d",(dim+1)%10u);
	sf_putint(mask,varname,ns);

	for (is=0; is < ns; is++) {
	    for (i=0; i < n123; i++) {
		kk[i] = msaa->mis[is][i]? 0.: 1.;
	    }
	    sf_intwrite (kk,n123,mask);
	}
    }

    msfind_pef (n123, dd, msaa, niter);

    for (i=0; i < dim; i++) {
	center[i] *= jump[0];
	a[i] *= jump[0];
    }

    onescale (0, msaa);
    print(dim, n, center, a, msaa->one);

    sf_putint(pef,"n1",nh);
    sf_putint(lag,"n1",nh);
    
    for (i=1; i < dim; i++) {
	sprintf(varname,"n%d",i+1);
	sf_putint(pef,varname,1);
	sf_putint(lag,varname,1);
    }

    sf_putint(lag,"n2",ns);
    sf_putints(lag,"jump",jump,ns);

    sf_intwrite(msaa->lag[0],nh*ns,lag);
    sf_floatwrite(msaa->flt,nh,pef);


    exit (0);
}

/* 	$Id$	 */
