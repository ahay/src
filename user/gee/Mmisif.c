/* Find MISSing Input values and Filter in 1-D. */
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

#include "misif.h"

int main(int argc, char* argv[])
{
    int n1, n2, nmiss, ia, na, lag, i2, im;
    float *xx, *aa;
    bool *mm;
    sf_file in, out, flt;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");
    flt = sf_output("filtout");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    n2 = sf_leftsize(in,1);

    xx = sf_floatalloc(n1);

    if (!sf_getint("nmiss",&nmiss)) nmiss=n1;
    /* number of iterations */

    if (!sf_getint("na",&na)) na=3;
    /* filter size */

    if (!sf_getint("lag",&lag)) lag=1;
    /* filter lag */

    aa = sf_floatalloc(na);
    mm = sf_boolalloc(n1+na);

    sf_putint(flt,"n1",na);
    sf_putint(flt,"lag",lag);

    for (i2=0; i2 < n2; i2++) { 
	for (ia=0; ia < na; ia++) {
	    aa[ia] = 0.;
	}
	aa[lag-1] = 1.;
	for (im=n1; im < n1+na; im++) {
	    mm[im] = false;
	}
	mm[n1+lag-1] = true;

	sf_floatread(xx,n1,in);

	for (im=0; im < n1; im++) {
	    mm[im] = (bool) (xx[im] != 0.);
	}

        misif1 (nmiss, na, n1, xx, aa, mm);

        sf_floatwrite (aa,na,flt);
        sf_floatwrite (xx,n1,out);
    }

    exit(0);
}

/* 	$Id$	 */
