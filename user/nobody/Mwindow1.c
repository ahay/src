/* Divide a dataset into smooth overlapping windows.

The input dimensions n1,... become n1,nw,...

If no taper parameter is present, running 
< window.rsf sfstack norm=no > input.rsf 
should reconstruct the original input.
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

#include "window1.h"

int main(int argc, char* argv[])
{
    int n2, i2, i1, n, nw, w, iw, i0;
    float h, *dat, *win, *dat2;
    bool taper, hastaper;
    sf_file in, out;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n)) sf_error("No n1= in input");
    n2 = sf_leftsize(in,1);

    if (!sf_getint ("nw",&nw)) sf_error("Need nw=");
    /* Number of windows */
    if (!sf_getint ("w",&w)) sf_error("Need w=");
    /* Window length */
    if (!sf_getfloat ("h",&h)) h = (nw*w - n)/(nw-1.);
    /* Window separation */

    hastaper = sf_getbool("taper",&taper);
    /* Window tapering */

    sf_putint(out,"n2",nw);
    sf_putint(out,"n3",n2);

    dat = sf_floatalloc (n);
    win = sf_floatalloc (w);
    dat2 = sf_floatalloc (n);
  
    window1_init (h,w-h);

    for (i2=0; i2 < n2; i2++) {
	sf_floatread (dat, n, in);
	for (iw=0; iw < nw; iw++) {
	    if (hastaper) {
		i0 = window1_apply(iw,w,dat,taper,taper,win);
	    } else {
		i0 = window1_apply(iw,w,dat,(iw > 0),(iw < nw-1),win);
	    }
	    for (i1=0; i1 < i0; i1++) {
		dat2[i1] = 0.;
	    }
	    for (i1=i0; i1 < i0+w; i1++) {
		dat2[i1] = win[i1-i0];
	    }
	    for (i1=i0+w; i1 < n; i1++) {
		dat2[i1] = 0.;
	    }
	    sf_floatwrite (dat2,n,out);
	}
    }

    exit (0);
}

/* 	$Id$	 */
