/* Cut an elliptic hole in data (for interpolation tests).
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

int main (int argc, char* argv[])
{
    int n1, n2, n3, i1, i2, i3;
    float *pp, x,y,u,v, *maskout;
    sf_file in, out, mask;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");
    mask = sf_output("maskout");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    n3 = sf_leftsize(in,2);

    pp = sf_floatalloc(n1);
    maskout = sf_floatalloc(n1);

    for (i3=0; i3 < n3; i3++) {
	for (i2=0; i2 < n2; i2++) { 	
	    sf_floatread (pp,n1,in);

	    for (i1=0; i1 < n1; i1++) { 
		x = ((float) i1)/n1 - 0.5;
		y = ((float) i2)/n2 - 0.3;
		u =  x+y;
		v = (x-y)/2.;
		if (u*u + v*v < 0.15 ) {
		    pp[i1] = 0.;
		    maskout[i1] = 0.;
		} else {
		    maskout[i1] = 1.;
		}
	    }
	    	    
	    sf_floatwrite (pp,n1,out);
	    sf_floatwrite (maskout,n1,mask);
	}
    }

    exit (0);
}

/* 	$Id$	 */
