/* Construct data from random lines */
/*
Copyright (C) 2000 The Board of Trustees of Stanford University
Copyright (C) 2010 University of Texas at Austin

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

int main(int argc, char* argv[])
{
    int i1, i2, n1,n2, iline, lines, mid1, mid2, i, seed;
    float **xx, rho, amp;
    sf_file out;

    sf_init(argc,argv);
    out = sf_output("out");
    sf_setformat(out,"native_float");

    if (!sf_getint("n1",&n1)) n1=10;
    if (!sf_getint("n2",&n2)) n2=10;
    /* dimensions */

    sf_putint(out,"n1",n1);
    sf_putint(out,"n2",n2);

    xx = sf_floatalloc2(n1,n2);

    if (!sf_getint("lines",&lines)) lines=3;
    /* number of lines */

    if (!sf_getint("seed",&seed)) seed=2000;
    /* random number seed */

    init_genrand((unsigned long) seed);

    for (i2=0; i2 < n2; i2++) {
	for (i1=0; i1 < n1; i1++) {
	    xx[i2][i1] = 0.;
	}
    }

    for (iline=0; iline < lines; iline++) {
	mid1 = n1 * genrand_real1();
	mid2 = n2 * genrand_real1();
	rho =  tanf( SF_PI/4 *(2*genrand_real1()-1.));
	amp =  2*genrand_real1()-1.;
	for (i=-n1; i <= n1; i++) {
	    i2 = i + mid1;
	    i1 = rho * i + mid2;
	    if(iline%2) {
		if (i2>=0 && i2 < n2 &&
		    i1>=0 && i1 < n1 ) xx[i2][i1] += amp;
	    } else {
		if (i2>=0 && i2 < n1 &&
		    i1>=0 && i1 < n2 ) xx[i1][i2] += amp;
	    }
	}
    }

    sf_floatwrite(xx[0],n1*n2,out);
    exit(0);
}

