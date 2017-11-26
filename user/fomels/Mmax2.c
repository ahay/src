/* Picking local maxima in 2-D */
/*
  Copyright (C) 2015 University of Texas at Austin
  
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

#include "max2.h"

int main(int argc, char* argv[])
{
    int n1, n2, n3, n12, i1, i2, i3, np, ip;
    float o1, o2, d1, d2, x, y;
    float **pick, **slice, **slice2;
    sf_file in, out;

    sf_init(argc, argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n1= in input");
    n3 = sf_leftsize(in,2);
    n12 = n1*n2;

    if (!sf_histfloat(in,"d1",&d1)) d1=1.;
    if (!sf_histfloat(in,"o1",&o1)) o1=0.;
    if (!sf_histfloat(in,"d2",&d2)) d2=1.;
    if (!sf_histfloat(in,"o2",&o2)) o2=0.;

    if (!sf_getint("np",&np)) np=n12;
    /* maximum number of picks */

    sf_putint(out,"n1",3);
    sf_putint(out,"n2",np);

    slice = sf_floatalloc2(n1,n2);
    pick = sf_floatalloc2(3,np);

    /* extend to avoid dealing with boundaries */
    slice2 = sf_floatalloc2(n1+2,n2+2);
    for (i2=0; i2 < n2+2; i2++) {
	slice2[i2][0] = slice2[i2][n1+1] = -SF_HUGE;
    }
    for (i1=0; i1 < n1+2; i1++) {
	slice2[0][i1] = slice2[n2+1][i1] = -SF_HUGE;
    }

    for (i3=0; i3 < n3; i3++) {
	sf_floatread(slice[0],n12,in);
	for (i2=0; i2 < n2; i2++) {
	    for (i1=0; i1 < n1; i1++) {
		slice2[i2+1][i1+1] = slice[i2][i1];
	    }
	}
	ip=0;
	for (i2=0; i2 < n2; i2++) {
	    for (i1=0; i1 < n1; i1++) {
		/* check for local maxima on a grid */
		if (not_max(slice2,i2+1,i1+1)) continue;

		/* now try to locate it more precisely */		
		pick[ip][2]=max2(slice2,i2,i1,&x,&y);
		pick[ip][0]=o1+(i1+x)*d1;
		pick[ip][1]=o2+(i2+y)*d2;

		ip++;
	    }
	}
	if (0==ip) {
	    pick[ip][0]=o1-d1;
	    pick[ip][1]=o2-d2;
	    pick[ip][2]=0.0f;
	    ip++;
	}
	for (i1=ip; i1 < np; i1++) {
	    pick[i1][0]=pick[ip-1][0];
	    pick[i1][1]=pick[ip-1][1];
	    pick[i1][2]=0.0f;
	}

	sf_floatwrite(pick[0],np*3,out);
    }

    exit(0);
}
