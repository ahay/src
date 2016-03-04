/* Recursive stacking by plane-wave construction. */
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
    int n1,n2,n3, i1,i2,i3, order;
    float eps, **u, **p, *trace;
    sf_file inp, out, dip;

    sf_init(argc,argv);
    inp = sf_input("in");
    dip = sf_input("dip");
    out = sf_output("out");

    if (!sf_histint(inp,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(inp,"n2",&n2)) sf_error("No n2= in input");
    n3 = sf_leftsize(inp,2);

    sf_putint(out,"n2",1);

    if (!sf_getbool("verb",&verb)) verb=false;
    if (!sf_getfloat("eps",&eps)) eps=0.01;
    /* regularization */
    
    if (!sf_getint("order",&order)) order=1;
    /* accuracy order */

    predict_init (n1, n2, eps*eps, order, 1, false);

    u = sf_floatalloc2(n1,n2);
    p = sf_floatalloc2(n1,n2);
    trace = sf_floatalloc(n1);

    for (i3=0; i3 < n3; i3++) {
	if (verb) sf_warning("cmp %d of %d;",i3+1,n3);

	sf_floatread(u[0],n1*n2,inp);
	sf_floatread(p[0],n1*n2,dip);

	/* load the last trace */
	for (i1=0; i1 < n1; i1++) {
	    trace[i1] = u[n2-1][i1];
	}
	for (i2=n2-2; i2 >= 0; i2--) {
	    predict_step(false,false,trace,p[i2]);
	    for (i1=0; i1 < n1; i1++) {
		trace[i1] += u[i2][i1];
	    }
	}

	/* normalize */
	for (i1=0; i1 < n1; i1++) {
	    trace[i1] /= n2;
	}

	sf_floatwrite(trace,n1,out);
    }
    if (verb) sf_warning(".");

    exit (0);
}
