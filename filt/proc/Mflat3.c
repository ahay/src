/* 3-D flattening (without picking). */
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
#include "predict.h"
#include "dijkstra.h"

int main (int argc, char *argv[])
{
    bool verb;
    int n1,n2,n3, n12, ref2, ref3, i2,i3,i1;
    float eps, ***u, ***p, ***q, **p2, **q2, ***v, pi, qi;
    sf_file in, out, idip, xdip;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");
    idip = sf_input("idip");
    xdip = sf_input("xdip");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    if (!sf_histint(in,"n3",&n3)) sf_error("No n3= in input");
    n12 = n1*n2*n3;

    if (!sf_getbool("verb",&verb)) verb=false;
    if (!sf_getfloat("eps",&eps)) eps=0.01;
    /* regularization */

    if (!sf_getint("ref2",&ref2)) ref2=0;
    if (!sf_getint("ref3",&ref3)) ref3=0;
    /* reference trace */

    u = sf_floatalloc3(n1,n2,n3);
    v = sf_floatalloc3(n1,n2,n3);
    p = sf_floatalloc3(n1,n2,n3);
    q = sf_floatalloc3(n1,n2,n3);

    sf_floatread(p[0][0],n12,idip);
    sf_floatread(q[0][0],n12,xdip);

    p2 = sf_floatalloc2(n2,n3);
    q2 = sf_floatalloc2(n2,n3);

    for (i3=0; i3 < n3; i3++) {
	for (i2=0; i2 < n2; i2++) {
	    pi = 0.;
	    qi = 0.;
	    for (i1=0; i1 < n1; i1++) {
		pi += p[i3][i2][i1]*p[i3][i2][i1]; 
		qi += q[i3][i2][i1]*q[i3][i2][i1];
	    }
	    p2[i3][i2] = pi;
	    q2[i3][i2] = qi;
	}
    }

    dijkstra_init(n2,n3);
    dijkstra(ref2,ref3,p2,q2);

    exit (0);
}

/* 	$Id: Mflat.c 743 2004-08-16 20:41:00Z fomels $	 */
