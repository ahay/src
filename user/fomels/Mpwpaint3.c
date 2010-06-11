/* 3-D painting by plane-wave construction. */
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
    int n1,n2,n3, n12, n23, ref2, ref3, i2,i3,i1, ud, lr, order;
    float eps, ***dat, ***p, ***q, **p2, **q2, *trace;
    sf_file dip, out, seed, cost;

    sf_init(argc,argv);
    dip = sf_input("in");
    out = sf_output("out");
    seed = sf_input("seed");
    cost = sf_input("cost");

    if (!sf_histint(dip,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(dip,"n2",&n2)) sf_error("No n2= in input");
    if (!sf_histint(dip,"n3",&n3)) sf_error("No n3= in input");
    n23 = n2*n3;
    n12 = n1*n23;

    if (!sf_getbool("verb",&verb)) verb=false;
    if (!sf_getfloat("eps",&eps)) eps=0.01;
    /* regularization */

    if (!sf_getint("ref2",&ref2)) ref2=0;
    if (!sf_getint("ref3",&ref3)) ref3=0;
    /* reference trace */

    sf_putint(out,"n4",1);

    p = sf_floatalloc3(n1,n2,n3);
    q = sf_floatalloc3(n1,n2,n3);
    dat = sf_floatalloc3(n1,n2,n3);

    sf_floatread(p[0][0],n12,dip);
    sf_floatread(q[0][0],n12,dip);

    p2 = sf_floatalloc2(n2,n3);
    q2 = sf_floatalloc2(n2,n3);

    sf_floatread(p2[0],n23,cost);
    sf_floatread(q2[0],n23,cost);

    dijkstra_init(n2,n3,p2,q2);
    dijkstra_source(ref2,ref3);
    sf_floatread(dat[ref3][ref2],n1,seed);
    if (verb) sf_warning("%d %d",ref2,ref3);

    if (!sf_getint("order",&order)) order=1;
    /* accuracy order */

    predict_init(n1,n2, eps*eps, order, 1, false);
    
    while (dijskstra_step(&i2,&i3,&ud,&lr)) {
	if (verb) sf_warning("%d %d",i2,i3);

	trace = dat[i3][i2];
	for (i1=0; i1 < n1; i1++) {
	    trace[i1] = dat[i3-lr][i2-ud][i1];
	} 
	if (ud > 0) {
	    predict_step(false,true,trace,p[i3][i2-ud]);
	} else if (ud < 0) {
	    predict_step(false,false,trace,p[i3][i2]);
	}

	if (lr > 0) {
	    predict_step(false,true,trace,q[i3-lr][i2]);
	} else if (lr < 0) {
	    predict_step(false,false,trace,q[i3][i2]);
	}
    }
	
    sf_floatwrite(dat[0][0],n12,out);

    exit (0);
}

/* 	$Id: Mflat.c 743 2004-08-16 20:41:00Z fomels $	 */
