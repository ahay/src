/* Wilson-Burg spectral factorization. */
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

#include "helix.h"
#include "wilson.h"
#include "compress.h"

int main(int argc, char* argv[])
{
    int ns, na, niter, maxlag, ia;
    float a0, s0, eps;
    filter ss, aa;
    char *lagfile;
    sf_file in, out, lag0, lag;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&ns)) sf_error("No n1= in input");
    ss = allocatehelix (ns);

    if (NULL == (lagfile = sf_histstring(in,"lag"))) {
	if (NULL == (lagfile = sf_getstring("lag"))) {
	    /* optional input file with filter lags */
	    for (ia=0; ia < ns; ia++) {
		ss->lag[ia]=ia+1;
	    }
	    lag0 = NULL;
	} else {
	    lag0 = sf_input("lag");
	}
    } else {
	lag0 = sf_input(lagfile);
    }

    if (NULL != lag0) {
	if (SF_INT != sf_gettype(lag0)) 
	    sf_error("Need int data in lag file '%s'",lagfile);

	sf_intread(ss->lag,ns,lag0);
    }
 
    maxlag = 0;
    for( ia=0; ia < ns; ia++) {
	if (ss->lag[ia] > maxlag) maxlag = ss->lag[ia];
    }

    if(!sf_histfloat(in,"a0",&s0)) s0 = 1.; 

    if (!sf_getint("niter",&niter)) niter=20;
    /* number of iterations */

    if (!sf_getfloat("eps",&eps)) eps=1.e-6;
    /* truncation tolerance */

    if (NULL == (lagfile = sf_getstring("lagin"))) {
	/* optional input file with output filter lags */
	if (!sf_getint("n1",&na)) na=maxlag;
	/* output filter length */

	aa = allocatehelix (na);

	for (ia=0; ia < na; ia++) {	    
	    aa->lag[ia]=ia+1;
	}
    } else {
	lag0 = sf_input("lagin");

	if (SF_INT != sf_gettype(lag0)) 
	    sf_error("Need int data in lag file '%s'",lagfile);

	if (!sf_histint(lag0,"n1",&na)) 
	    sf_error("No n1= in lag file '%s'",lagfile);

	aa = allocatehelix (na);

	sf_intread(aa->lag,na,lag0);
    }
    for (ia=0; ia < na; ia++) {	    
	aa->flt[ia]=0.;
    }

    sf_floatread(ss->flt,ns,in);

    wilson_init( maxlag*10);
    a0 = wilson_factor(niter, 2.*s0, ss, aa, true, 1.e-6);

    aa = compress(aa,eps);
    na = aa->nh;

    lag = sf_output("lagout");
    sf_putint(lag,"n1",na);
    sf_settype(lag,SF_INT);
    sf_fileflush(lag,lag0);

    sf_intwrite(aa->lag,na,lag);
    sf_fileclose(lag);

    if (NULL != (lagfile = sf_getstring("lagout"))) 
	sf_putstring(out,"lag",lagfile);

    sf_putint(out,"n1",na);
    sf_putfloat(out,"a0",a0);

    for( ia=0; ia < na; ia++) {
	aa->flt[ia] *= a0;
    }
    sf_floatwrite(aa->flt,na,out);

    exit (0);
}

/* 	$Id$	 */
