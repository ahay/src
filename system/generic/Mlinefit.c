/* Fit a line to a set of points in 2-D.
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
#include "spline3.h"

int main(int argc, char* argv[])
{
    int id, nd, n1, i1, n2, i2, two;
    float **table=NULL, *trace=NULL;
    float x, y, o1, d1, a, b, sx, sx2, sxy, sy, det;
    sf_file in=NULL, out=NULL, pattern=NULL;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");

    if (!sf_histint(in,"n1",&two) || two != 2) sf_error("Need n1=2 in input");
    if (!sf_histint(in,"n2",&nd)) sf_error ("Need n2= in input");
    sf_putint(out,"n2",1);
    n2 = sf_leftsize(in,2);

    if (NULL != sf_getstring("pattern")) {
	pattern = sf_input("pattern");
    } else {
	pattern = NULL;
    }

    if (!sf_getint("n1",&n1) &&
	(NULL== pattern ||
	 !sf_histint(pattern,"n1",&n1))) sf_error("Need n1=");
    /* Output grid size */
    if (!sf_getfloat("d1",&d1) &&
	(NULL== pattern ||
	 !sf_histfloat(pattern,"d1",&d1))) sf_error("Need d1=");
    /* Output sampling */
    if (!sf_getfloat("o1",&o1) &&
	(NULL== pattern ||
	 !sf_histfloat(pattern,"o1",&o1))) sf_error("Need o1=");
    /* Output origin */

    sf_putint(out,"n1",n1);
    sf_putfloat(out,"o1",o1);
    sf_putfloat(out,"d1",d1);

    trace = sf_floatalloc(n1);
    table = sf_floatalloc2(2,nd);

    for (i2=0; i2 < n2; i2++) {
	sf_floatread(table[0],2*nd,in);

	sx = sx2 = sxy = sy = 0.;
	for (id=0; id < nd; id++) {
	    x = table[id][0];
	    y = table[id][1];
	    sx += x;
	    sx2 += x*x;
	    sxy += x*y;
	    sy += y;
	}
	det = nd*sx2-sx*sx;

	if (0.==det) sf_error("zero determinant at trace %d",i2);

	a = (nd*sxy-sx*sy)/det;
	b = (sy*sx2-sx*sxy)/det;

	for (i1=0; i1 < n1; i1++) {
	    x = o1 + i1*d1;
	    trace[i1] = a*x+b;
	}

	sf_floatwrite(trace,n1,out);
    }

    exit(0);
}
