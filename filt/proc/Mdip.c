/* 3-D dip estimation by plane wave destruction.
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

#include "dip3.h"

int main (int argc, char *argv[])
{
    int n1,n2,n3, n123, niter, nw, nj1, nj2, i, rect[3], liter;
    float p0, q0, ***u, ***p;
    bool verb, sign;
    sf_file in, out;

    sf_init(argc,argv);
    in = sf_input ("in");
    out = sf_output ("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float type");

    if (!sf_histint(in,"n1",&n1)) sf_error("Need n1= in input");
    if (!sf_histint(in,"n2",&n2)) n2=1;
    if (!sf_histint(in,"n3",&n3)) n3=1;
    n123 = n1*n2*n3;

    /* two dips output in 3-D */
    if (n3 > 1) sf_putint(out,"n4",2); 

    if (!sf_getint("niter",&niter)) niter=5;
    /* number of iterations */
    if (!sf_getint("liter",&liter)) liter=20;
    /* number of linear iterations */

    if (!sf_getint("rect1",&rect[0])) rect[0]=1;
    if (!sf_getint("rect2",&rect[1])) rect[1]=1;
    if (!sf_getint("rect3",&rect[2])) rect[2]=1;
    /* dip smoothness */

    if (!sf_getfloat("p0",&p0)) p0=0.;
    /* initial in-line dip */
    if (!sf_getfloat("q0",&q0)) q0=0.;
    /* initial cross-line dip */

    if (!sf_getint("order",&nw)) nw=1;
    /* [1,2,3] accuracy order */
    if (nw < 1 || nw > 3) 
	sf_error ("Unsupported nw=%d, choose between 1 and 3",nw);
    if (!sf_getint("nj1",&nj1)) nj1=1;
    /* in-line antialiasing */
    if (!sf_getint("nj2",&nj2)) nj2=1;
    /* cross-line antialiasing */

    if (!sf_getbool("verb",&verb)) verb = false;
    /* verbosity flag */
    if (!sf_getbool("sign",&sign)) sign = false;
    /* if y, keep dip sign constant */
    
    /* initialize dip estimation */
    dip3_init(n1, n2, n3, rect, liter, sign);

    u = sf_floatalloc3(n1,n2,n3);
    p = sf_floatalloc3(n1,n2,n3);

    /* read data */
    sf_floatread(u[0][0],n123,in);

    /* initialize t-x dip */
    for(i=0; i < n123; i++) {
	p[0][0][i] = p0;
    }
  
    /* estimate t-x dip */
    dip3(1, niter, nw, nj1, verb, u, p);

    /* write t-x dip */
    sf_floatwrite(p[0][0],n123,out);

    if (n3 > 1) { /* if 3-D input */
	/* initialize t-y dip */
	for(i=0; i < n123; i++) {
	    p[0][0][i] = q0;
	}
  
	/* estimate t-y dip */
	dip3(2, niter, nw, nj2, verb, u, p);

	/* write t-y dip */
	sf_floatwrite(p[0][0],n123,out);
    }
    
    sf_close();
    exit (0);
}

/* 	$Id: Mdip.c,v 1.7 2004/06/25 18:08:42 fomels Exp $	 */
