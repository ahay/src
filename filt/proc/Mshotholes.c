/* Remove random shot gathers from a 2-D dataset.
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

#include <float.h>

#include <rsf.h>

#include "randn.h"

int main(int argc, char* argv[])
{
    int n1, n2, n3, i1, i2, i3, is, **known;
    float *trace, *zero, *chance, perc, rand, sum;
    sf_file in, out, mask;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");
    mask = sf_output("mask");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    if (!sf_histint(in,"n3",&n3)) sf_error("No n3= in input");

    sf_putint(mask,"n1",n2);
    sf_putint(mask,"n2",n3);
    sf_putint(mask,"n3",1);
    sf_setformat(mask,"native_int");

    if (!sf_getfloat("perc",&perc)) perc=0.75;
    /* how many shots to remove */

    trace = sf_floatalloc(n1);
    zero = sf_floatalloc(n1);
    known = sf_intalloc2(n2,n3);
    chance = sf_floatalloc(n2+n3);

    for (i1=0; i1 < n1; i1++) {
	zero[i1] = 0.;
    }

    random0 (n2+n3,chance);

    for (i3=0; i3 < n3; i3++) { /* half-offset */
	for (i2=0; i2 < n2; i2++) { /* midpoint */
	    is = i2 - i3 + n3-1; /* shot */
	    rand = chance[is];

	    sf_floatread (trace,n1,in);
	    sum = 0.;
	    for (i1=0; i1 < n1; i1++) {
		sum += trace[i1]*trace[i1];
	    }

	    if (rand > perc && sum > FLT_EPSILON) {
		sf_floatwrite (trace,n1,out);
		known[i3][i2] = 1;
	    } else {
		sf_floatwrite (zero,n1,out);
		known[i3][i2] = 0;
	    }
	}
    }
    sf_intwrite (known[0],n2*n3,mask);

    exit(0);
}

/* 	$Id: Mshotholes.c,v 1.6 2004/07/02 11:54:48 fomels Exp $	 */
