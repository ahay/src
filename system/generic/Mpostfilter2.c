/* Convert B-spline coefficients to data in 2-D. */
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

int main(int argc, char* argv[])
{
    int n1, n2, n3, n12, i3, i2, i1, nw;
    float *slice=NULL, *slice2=NULL;
    bool vert, horz;
    sf_file in=NULL, out=NULL;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    n12 = n1*n2;
    n3 = sf_leftsize(in,2);

    if (!sf_getint("nw",&nw)) sf_error("Need nw=");
    /* filter size */

    if (!sf_getbool("vert",&vert)) vert=true;
    /* include filter on the first axis */
    if (!sf_getbool("horz",&horz)) horz=true;
    /* include filter on the second axis */

    slice = sf_floatalloc(n12);
    slice2 = sf_floatalloc(n12);

    for (i3=0; i3 < n3; i3++) {
	sf_floatread(slice,n12,in);

	if (nw > 2) { 
	    if (vert) {
		for (i2=0; i2 < n2; i2++) {
		    sf_spline_post (nw, i2*n1, 1, n1, slice, slice2);
		}
	    } else {
		for (i1=0; i1 < n12; i1++) {
		    slice2[i1] = slice[i1];
		}
	    }

	    if (horz) {
		for (i1=0; i1 < n1; i1++) {
		    sf_spline_post (nw, i1, n1, n2, slice2, slice);
		}
	    } else {
		for (i1=0; i1 < n12; i1++) {
		    slice[i1] = slice2[i1];
		}
	    }
	}

	sf_floatwrite(slice,n12,out);
    }

    exit(0);
}
