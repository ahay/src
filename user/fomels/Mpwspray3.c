/* Plane-wave spray. */
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

int main (int argc, char *argv[])
{
    bool verb;
    int n1,n2,n3, i1,i2,i3, is1, ns1, is2, ns2, ip, np1, np2, np;
    float eps, ****u, ***p1, ***p2, **c1, **c2, *trace;
    sf_file inp, out, dip, cost;

    sf_init(argc,argv);
    inp = sf_input("in");
    dip = sf_input("dip");
    cost = sf_input("cost");
    out = sf_output("out");

    if (!sf_histint(dip,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(dip,"n2",&n2)) sf_error("No n2= in input");
    if (!sf_histint(dip,"n3",&n3)) sf_error("No n3= in input");

    if (!sf_getbool("verb",&verb)) verb=false;
    /* verbosity */
    if (!sf_getfloat("eps",&eps)) eps=0.01;
    /* regularization */
    
    if (!sf_getint("ns1",&ns1)) sf_error("Need ns1=");
    if (!sf_getint("ns2",&ns2)) sf_error("Need ns2=");
    /* spray radius */
    np1 = 2*ns1+1;
    np2 = 2*ns2+1;
    np = np1*np2;

    sf_putint(out,"n2",np);
    sf_shiftdim(inp, out, 2);

    c1 = sf_floatalloc2(n2,n3);
    c2 = sf_floatalloc2(n2,n3);

    sf_floatread(c1[0],n2*n3,cost);
    sf_floatread(c2[0],n2*n3,cost);

    u = sf_floatalloc4(n1,np,n2,n3);
    for (i3=0; i3 < n3; i3++) {
	for (i2=0; i2 < n2; i2++) {
	    for (ip=0; ip < np; ip++) {
		for (i1=0; i1 < n1; i1++) {
		    u[i3][i2][ip][i1] = 0.;
		}
	    }
	}
    }

    p1 = sf_floatalloc3(n1,n2,n3);
    p2 = sf_floatalloc3(n1,n2,n3);
    trace = sf_floatalloc(n1);

    sf_floatread(p1[0][0],n1*n2*n3,dip);
    sf_floatread(p2[0][0],n1*n2*n3,dip);

    for (i3=0; i3 < n3; i3++) { 
	for (i2=0; i2 < n2; i2++) { 
	    sf_floatread(u[i3][i2][ns2*np1+ns1],n1,inp);


	    dijkstra_init(n2,n3,p2,q2);
	    dijkstra_source(i2,i3);

	    while (dijskstra_step(&k2,&k3,&ud,&lr)) {
		
		/* predict k3,k2 from k3-lr,k2-ud */

		trace = dat[k3][k2][(k3-i3+
;
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

	    /* predict forward */
	    for (is=0; is < ns; is++) {
		ip = i2-is-1;
		if (ip < 0) break;
		predict_step(false,false,trace,p[ip]);
		for (i1=0; i1 < n1; i1++) {
		    u[ip][ns-is-1][i1] = trace[i1];
		}
	    }

	}
	sf_floatwrite(u[0][0],n1*ns2*n2,out);
    }	    

    exit (0);
}

/* 	$Id: Mflat.c 1131 2005-04-20 18:19:10Z fomels $	 */
