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

#include "spline.h"

int main(int argc, char* argv[])
{
    int n1, n2, n3, n12, i3, i2, i1, nw;
    float *slice, *slice2;
    sf_file in, out;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    n12 = n1*n2;
    n3 = sf_leftsize(in,2);

    if (!sf_getint("nw",&nw)) sf_error("Need nw=");
    /* filter size */

    slice = sf_floatalloc(n12);
    slice2 = sf_floatalloc(n12);

    for (i3=0; i3 < n3; i3++) {
	sf_floatread(slice,n12,in);

	if (nw > 2) { 
	    for (i2=0; i2 < n2; i2++) {
		spline_post (nw, i2*n1, 1, n1, slice, slice2);
	    }
	    for (i1=0; i1 < n1; i1++) {
		spline_post (nw, i1, n1, n2, slice2, slice);
	    }
	}

	sf_floatwrite(slice,n12,out);
    }
    
    exit(0);
}

