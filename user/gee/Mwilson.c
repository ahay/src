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

#include "wilson.h"
#include "wilson2.h"
#include "compress.h"

int main(int argc, char* argv[])
{
    int ns, na, niter, maxlag, ia;
    float a0, s0, eps;
    sf_filter ss, aa;
    bool verb, stable;
    char *file;
    sf_file in, out, lag0, lag, flt0;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&ns)) sf_error("No n1= in input");
    ss = sf_allocatehelix (ns);

    if (NULL == (file = sf_histstring(in,"lag"))) {
	if (NULL == (file = sf_getstring("lag"))) {
	    /* optional input file with filter lags */
	    for (ia=0; ia < ns; ia++) {
		ss->lag[ia]=ia+1;
	    }
	    lag0 = NULL;
	} else {
	    lag0 = sf_input("lag");
	}
    } else {
	lag0 = sf_input(file);
    }

    if (NULL != lag0) {
	if (SF_INT != sf_gettype(lag0)) 
	    sf_error("Need int data in lag file '%s'",file);

	sf_intread(ss->lag,ns,lag0);
    }
 
    if (!sf_getint("maxlag",&maxlag)) {
	/* maximum lag */
	maxlag = 0;
	for( ia=0; ia < ns; ia++) {
	    if (ss->lag[ia] > maxlag) maxlag = ss->lag[ia];
	}
    }

    if(!sf_histfloat(in,"a0",&s0)) s0 = 1.; 

    if (!sf_getint("niter",&niter)) niter=20;
    /* number of iterations */

    if (!sf_getfloat("eps",&eps)) eps=1.e-6;
    /* truncation tolerance */

    if (NULL == (file = sf_getstring("lagin"))) {
	/* optional input file with output filter lags */
	if (!sf_getint("n1",&na)) na=maxlag;
	/* output filter length */

	aa = sf_allocatehelix (na);

	for (ia=0; ia < na; ia++) {	    
	    aa->lag[ia]=ia+1;
	}
    } else {
	lag0 = sf_input("lagin");

	if (SF_INT != sf_gettype(lag0)) 
	    sf_error("Need int data in lag file '%s'",file);

	if (!sf_histint(lag0,"n1",&na)) 
	    sf_error("No n1= in lag file '%s'",file);

	aa = sf_allocatehelix (na);

	sf_intread(aa->lag,na,lag0);
    }

    if (!sf_getfloat("a0",&a0)) a0=1.;

    if (NULL == (file = sf_getstring("filtin"))) {
	for (ia=0; ia < na; ia++) {	    
	    aa->flt[ia]=0.;
	}
    } else {
	flt0 = sf_input(file);
	sf_floatread(aa->flt,na,flt0);
	sf_fileclose(flt0);
    }

    if (!sf_getbool("verb",&verb)) verb=true;
    /* verbosity flag */
    if (!sf_getbool("stable",&stable)) stable=false;
    /* stability flag */

    sf_floatread(ss->flt,ns,in);

    if (stable) {
	wilson2_init( maxlag*10);
	a0 = wilson2_factor(niter, 2.*s0, ss, a0, aa, verb, 1.e-6);
    } else {
	wilson_init( maxlag*10);
	a0 = wilson_factor(niter, 2.*s0, ss, aa, verb, 1.e-6);
    }

    aa = compress(aa,eps);
    na = aa->nh;

    lag = sf_output("lagout");
    sf_putint(lag,"n1",na);
    sf_settype(lag,SF_INT);
    sf_fileflush(lag,lag0);

    sf_intwrite(aa->lag,na,lag);
    sf_fileclose(lag);

    if (NULL != (file = sf_getstring("lagout"))) 
	sf_putstring(out,"lag",file);

    sf_putint(out,"n1",na);
    sf_putfloat(out,"a0",a0);

    for( ia=0; ia < na; ia++) {
	aa->flt[ia] *= a0;
    }
    sf_floatwrite(aa->flt,na,out);


    exit (0);
}

/* 	$Id: Mwilson.c 7107 2011-04-10 02:04:14Z ivlad $	 */
