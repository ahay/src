/* Autocorrelation for helix filters. */
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

#include "autocorr.h"

int main(int argc, char* argv[])
{
    int i, na, ns;
    float s0, a0;
    char* lagfile;
    sf_filter ss, aa;
    sf_file in, out, lag0, lag;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&na)) sf_error("No n1= in input");
    aa = sf_allocatehelix (na);

    if (NULL == (lagfile = sf_histstring(in,"lag"))) {
	if (NULL == (lagfile = sf_getstring("lag"))) {
	    /* optional input file with filter lags */
	    for (i=0; i < na; i++) {
		aa->lag[i]=i+1;
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

	sf_intread(aa->lag,na,lag0);
    }

    if (!sf_histfloat(in,"a0",&a0)) a0=1.;
    sf_floatread (aa->flt,na,in);

    ss = autocorr (aa, a0, &s0, 1.e-6);
    ns = ss->nh;

    sf_putfloat (out,"a0",s0*0.5);
    sf_putint (out,"n1",ns);

    lag = sf_output("lagout");
    sf_putint(lag,"n1",ns);
    sf_settype(lag,SF_INT);
    sf_fileflush(lag,lag0);

    sf_intwrite(ss->lag,ns,lag);
    sf_fileclose(lag);

    if (NULL != (lagfile = sf_getstring("lagout"))) 
	sf_putstring(out,"lag",lagfile);

    sf_floatwrite (ss->flt,ns,out);


    exit (0);
}

/* 	$Id: Mautocorr.c 7107 2011-04-10 02:04:14Z ivlad $	 */
