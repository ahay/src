/* Cosine taper around the borders (2-D).
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

#include <rsf.h>

int main(int argc, char* argv[])
{
    int n1, n2, nw1, nw2, i1, i2, iw;
    float *trace, *w1=NULL, *w2=NULL, wi;
    sf_file in, out;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    n2 = sf_leftsize(in,1);

    if (!sf_getint("nw1",&nw1)) nw1=0;
    /* tapering length on first axis */
    if (nw1 > n1) nw1=n1;

    if (!sf_getint("nw2",&nw2)) nw2=0;
    /* tapering length on second axis */
    if (nw2 > n2) nw2=n2;

    trace = sf_floatalloc(n1);
    if (nw1 > 0) w1 = sf_floatalloc(nw1);
    if (nw2 > 0) w2 = sf_floatalloc(nw2);

    for (iw=0; iw < nw1; iw++) {
	wi = sinf(0.5*SF_PI*(iw+1.)/(nw1+1.));
	w1[iw] = wi*wi;
    }

    for (iw=0; iw < nw2; iw++) {
	wi = sinf(0.5*SF_PI*(iw+1.)/(nw2+1.));
	w2[iw] = wi*wi;
    }

    for (i2=0; i2 < n2; i2++) {
	sf_floatread(trace,n1,in);
	for (iw=0; iw < nw1; iw++) {
	    wi = w1[iw];
	    trace[iw]      *= wi;
	    trace[n1-1-iw] *= wi;
	}
	if (i2 < nw2) {
	    wi = w2[i2];
	    for (i1=0; i1 < n1; i1++) {
		trace[i1] *= wi;
	    }
	}
	if (n2-1-i2 < nw2) {
	    wi = w2[n2-1-i2];
	    for (i1=0; i1 < n1; i1++) {
		trace[i1] *= wi;
	    }
	}
	sf_floatwrite(trace,n1,out);
    }

    sf_close();
    exit(0);
}

/* 	$Id: Mcostaper.c,v 1.4 2004/06/25 08:41:19 fomels Exp $	 */
