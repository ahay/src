/* Canny-like edge detector. */
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

#include "edge.h"

int main(int argc, char* argv[])
{
    int n1, n2, n3, i1, i2, i3, j1, j2, k1, k2;
    float **pp, **ww, **w1, **w2, g1, g2;
    sf_file in, out;

    sf_init(argc,argv);

    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    n3 = sf_leftsize(in,2);

    pp = sf_floatalloc2(n1,n2);
    w1 = sf_floatalloc2(n1,n2);
    w2 = sf_floatalloc2(n1,n2);
    ww = sf_floatalloc2(n1,n2);

    for (i3=0; i3 < n3; i3++) {
	sf_floatread(pp[0],n1*n2,in);
	/* gradient computation */
	grad3(n1,n2,pp,w1,w2);
	for (i2=0; i2 < n2; i2++) {
	    for (i1=0; i1 < n1; i1++) {
		/* gradient norm */
		g1 = w1[i2][i1];
		g2 = w2[i2][i1];
		ww[i2][i1] = g1*g1+g2*g2;
	    }
	}
	/* edge thinning */
	for (i2=0; i2 < n2; i2++) {
	    for (i1=0; i1 < n1; i1++) {
		g1 = w1[i2][i1];
		g2 = w2[i2][i1];
		if (fabsf(g1) > fabsf(g2)) {
		    j1=1;
		    if (g2/g1 > 0.5) { 
			j2=1;
		    } else if (g2/g1 < - 0.5) {
			j2=-1; 
		    } else {
			j2=0;
		    }
		} else if (fabsf(g2) > fabsf(g1)) {
		    j2=1;
		    if (g1/g2 > 0.5) { 
			j1=1;
		    } else if (g1/g2 < - 0.5) {
			j1=-1; 
		    } else {
			j1=0;
		    }
		} else {
		    j1=0;
		    j2=0;
		}
		k1 = i1+j1; if (j1 && (k1 < 0 || k1 >= n1)) k1=i1;
		k2 = i2+j2; if (j2 && (k2 < 0 || k2 >= n2)) k2=i2;
		if (ww[i2][i1] <= ww[k2][k1]) {
		    pp[i2][i1] = 0.;
		    continue;
		} 
		k1 = i1-j1; if (k1 < 0 || k1 >= n1) k1=i1;
		k2 = i2-j2; if (k2 < 0 || k2 >= n2) k2=i2;
		if (ww[i2][i1] <= ww[k2][k1]) {
		    pp[i2][i1] = 0.;
		    continue;
		} 
		pp[i2][i1] = ww[i2][i1];
	    }
	}
	/* edge selection */

	sf_floatwrite(pp[0],n1*n2,out);
    }

    exit(0);
}
