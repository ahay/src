/* Divide a dataset into 2-D smooth overlapping windows.
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

#include "window2.h"

int main(int argc, char* argv[])
{
    int n[2], nw[2], w[2], i0[2], i[2], i1, i2, j;
    float h[2], **dat, **win, **dat2;
    sf_file in, out;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",n))   sf_error("No n1= in input");
    if (!sf_histint(in,"n2",n+1)) sf_error("No n2= in input");

    if (!sf_getints("nw",nw,2)) sf_error("Need nw=");
    /* number of windows */
    if (!sf_getints("w",w,2)) sf_error("Need w=");
    /* window size */
    if (!sf_getfloats ("h",h,2)) {
	for (j=0; j < 2; j++) {
	    h[j] = (nw[j]*w[j] - n[j])/(nw[j]-1.);
	}
    }
    sf_putint(out,"n3",nw[0]*nw[1]);

    dat = sf_floatalloc2(n[0],n[1]);
    win = sf_floatalloc2(w[0],w[1]);
    dat2 = sf_floatalloc2(n[0],n[1]);
  
    window2_init (w,nw,n,h);

    sf_floatread (dat[0],n[0]*n[1],in);

    for (i[1]=0; i[1] < nw[1]; i[1]++) {
	for (i[0]=0; i[0] < nw[0]; i[0]++) {
	    window2_apply(i,dat,
			  (i[1] > 0),(i[1] < nw[1]-1),
			  (i[0] > 0),(i[0] < nw[0]-1),i0,win);
	    for (i2=0; i2 < n[1]; i2++) {
		for (i1=0; i1 < n[0]; i1++) {
		    if (i1 >= i0[0] && i1 < i0[0]+w[0] &&
			i2 >= i0[1] && i2 < i0[1]+w[1]) {
			dat2[i2][i1] = win[i2-i0[1]][i1-i0[0]];
		    } else {
			dat2[i2][i1] = 0.;
		    }
		}
	    }
	    sf_floatwrite (dat2[0],n[0]*n[1],out);
	}
    }

    exit (0);
}

/* 	$Id$	 */
