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
#include <math.h>

#include <rsf.h>

#include "predict.h"
#include "dijkstra.h"

int main (int argc, char *argv[])
{
    bool verb;
    int n1,n2,n3, n12, ref2, ref3, i2,i3,i1, ud, lr, j3, j2;
    float eps, *trace, ***p, ***q, **p2, **q2, pi, qi;
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
	    p2[i3][i2] = sqrtf(pi/n1);
	    q2[i3][i2] = sqrtf(qi/n1);
	}
    }

    dijkstra_init(n2,n3);
    dijkstra(ref2,ref3,p2,q2);

    trace = sf_floatalloc(n1);
    predict_init(n1,n2,eps,1);

    for (i3=0; i3 < n3; i3++) {
	for (i2=0; i2 < n2; i2++) {
	    sf_floatread(trace,n1,in);
	    dijkstra_start(i2,i3);
	    j2 = i2;
	    j3 = i3;
	    if (verb) sf_warning("start with %d,%d",j2,j3);
	    while (dijkstra_next(&ud,&lr)) {
		if (0==lr) {
		    if (ud > 0) {
			j2 -= ud;
			predict_step(false,false,trace,p[j3][j2]);
		    } else {
			predict_step(false,true,trace,p[j3][j2]);
			j2 -= ud;
		    }
		} else if (0==ud) {
		    if (lr > 0) {
			j3 -= lr;
			predict_step(false,false,trace,q[j3][j2]);
		    } else {
			predict_step(false,true,trace,q[j3][j2]);
			j3 -= lr;
		    }
		}

		if (verb) sf_warning("then %d,%d",j2,j3);
	    }
	    sf_floatwrite(trace,n1,out);
	}
    }

    exit (0);
}

/* 	$Id: Mflat.c 743 2004-08-16 20:41:00Z fomels $	 */
