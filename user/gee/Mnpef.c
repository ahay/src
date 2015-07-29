/* Estimate Non-stationary PEF in N dimensions.
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

#include <math.h>
#include <float.h>

#include <rsf.h> 

#include "nhelix.h"
#include "createnhelix.h"
#include "nmisinput.h"
#include "npef.h"
#include "random.h"
 
int main(int argc, char* argv[])
{
    int n[SF_MAX_DIM], a[SF_MAX_DIM], center[SF_MAX_DIM], gap[SF_MAX_DIM];
    int *pch, *nh, dim, n123, nf, i, niter, nbf, nbp, id, ip, ig, np;
    int *kk, *pp;
    float *dd, eps, dabs, di;
    nfilter aa, bb;
    char varname[6], *lagfile;
    sf_file in, flt, lag, mask, patch, reg;

    sf_init(argc,argv);
    in = sf_input("in");
    flt = sf_output("out");

    dim = sf_filedims(in,n);
    
    if (NULL == (lagfile = sf_getstring("lag"))) sf_error("Need lag=");
    /* output file for filter lags */

    lag = sf_output(lagfile);
    sf_settype(lag,SF_INT);

    sf_putstring(flt,"lag",lagfile);

    sf_putints(lag,"n",n,dim);

    if (!sf_getints("a",a,dim)) sf_error("Need a=");

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
	    kk[i] = 1;
	}
    }

    sf_floatread(dd,n123,in);

    dabs = fabsf(dd[0]);
    for (i=1; i < n123; i++) {
	di = fabsf(dd[i]);
	if (di > dabs) dabs=di;
    }

    random_init(2004);
    for (i=0; i < n123; i++) {
	dd[i] = dd[i]/dabs+ 100.*FLT_EPSILON*(random0()-0.5);;
    }

    pp = sf_intalloc(n123);
    if (NULL != sf_getstring("pch")) {
	patch = sf_input("pch");
	if (SF_INT != sf_gettype(patch)) sf_error("Need int pch");

	sf_intread(pp,n123,patch);

	np = pp[0];
	for (i=1; i < n123; i++) {
	    if (pp[i] > np) np = pp[i];
	}

	sf_fileclose(patch);
    } else {
	np = n123;
	for (i=0; i < n123; i++) {
	    pp[i] = i;
	}
    }

    aa = createnhelix(dim, n, center, gap, a, pp);
    free (pp);

    nf = aa->hlx[0]->nh;
    nfind_mask(n123, kk, aa);

    if(!sf_getint("niter",&niter)) niter=100;
    /* number of iterations */
    if (!sf_getfloat("epsilon",&eps)) eps=0.01;
    /* regularization parameter */

    sf_putint(flt,"n1",nf);
    sf_putint(flt,"n2",np);

    sf_putint(lag,"n1",nf);
    sf_putint(lag,"n2",np);

    for (i=2; i < dim; i++) {
	sprintf(varname,"n%d",i+1);
	sf_putint(flt,varname,1);
	sf_putint(lag,varname,1);
    }

    for (ip=0; ip < np; ip++) {
	sf_intwrite(aa->hlx[ip]->lag,nf,lag);
    }
    sf_fileclose(lag);

    if (NULL != sf_getstring("maskout")) {
	/* optional output mask file */
	mask = sf_output("maskout");

	for (i=0; i < n123; i++) {
	    kk[i] = aa->mis[i]? 0.: 1.;
	}
	
	sf_settype(mask,SF_INT);
	sf_intwrite (kk,n123,mask);
    }
    
    reg = sf_input("filt");
    if (!sf_histint(reg,"n1",&nbf)) sf_error("No n1= in filt");
    if (!sf_histint(reg,"n2",&nbp)) sf_error("No n2= in filt");
    
    if (NULL != sf_getstring("filt_pch")) {
	patch = sf_input("filt_pch");
	if (SF_INT != sf_gettype(patch)) sf_error("Need int filt_pch");


	pp = sf_intalloc(np);
	sf_intread(pp,np,patch);
    } else {
	if (nbp != np) sf_error ("Wrong filter size: %d != %d",nbp,np);
	pp = NULL;
    }

    pch = sf_intalloc(nf*np);
    nh = sf_intalloc(nbp);

    for (i=0; i < nbp; i++) {
	nh[i] = nbf;
    }
    
    for (id=ig=0; ig < nf; ig++) {
	for (ip=0; ip < np; ip++, id++) {
	    pch[id] = (NULL != pp)? pp[ip]: ip;
	}
    }

    bb = nallocate (nbp, nf*np, nh, pch);

    if (NULL == (lagfile = sf_getstring("filt_lag")) &&
	NULL == (lagfile = sf_histstring(reg,"lag"))) 
	sf_error("Need filt_lag=");
    /* input file for double-helix filter lags */

    lag = sf_input(lagfile);
    if (SF_INT != sf_gettype(lag)) sf_error("Need int filt_lag");

    for (ip=0; ip < nbp; ip++) {
	sf_intread (kk,nbf,lag);
	for (i=0; i < nbf; i++) {
	    bb->hlx[ip]->lag[i] = kk[i]*nf;
	}
    }

    for (ip=0; ip < nbp; ip++) {
	sf_floatread (bb->hlx[ip]->flt,nbf,reg);
    }

    nfind_pef (n123, dd, aa, bb, niter, eps, nf);

    for (ip=0; ip < np; ip++) {
	sf_floatwrite (aa->hlx[ip]->flt,nf,flt);
    }


    exit(0);
}

/* 	$Id: Mnpef.c 11207 2013-10-28 16:56:21Z sfomel $	 */
