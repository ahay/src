/* Missing data interpolation in 2-D by Laplacian regularization. */
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
#include "lapfill.h"

int main(int argc, char* argv[])
{
    int n1,n2, n12, n3, i1, i3, niter;
    float *map, *msk;
    bool *known, grad, verb;
    sf_file in, out, mask;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    n3 = sf_leftsize(in,2);
    n12 = n1*n2;

    if (!sf_getint("niter",&niter)) niter=200;
    /* number of iterations */
    if (!sf_getbool("grad",&grad)) grad=false;
    /* if y, use gradient instead of laplacian */
    if (!sf_getbool("verb",&verb)) verb=false;
    /* verbosity flag */

    map = sf_floatalloc(n12);
    known = sf_boolalloc(n12);

    if (NULL != sf_getstring("mask")) {
	/* optional mask file with zeroes for missing data locations */
	mask = sf_input("mask");
	msk =  sf_floatalloc(n12);
    } else {
	mask = NULL;
	msk = map;
    }

    lapfill_init (n1,n2,grad,verb);

    for (i3=0; i3 < n3; i3++) {
	sf_floatread(map,n12,in);
	
	if (NULL != mask) sf_floatread(msk,n12,mask);

	for (i1=0; i1 < n12; i1++) {
	    known[i1] = (bool) (msk[i1] != 0.);
	}

	lapfill(niter,map,known);

	sf_floatwrite(map,n12,out);
    }

    exit(0);
}
