/* Plane-wave smoothing by box cascade. */
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
    int n1,n2,n3, n12, i1,i2,i3, order, i, n;
    float **u, **p, **next, *left, *right, *trace;
    float eps, r;
    sf_file inp, out, dip;

    sf_init(argc,argv);
    inp = sf_input("in");
    dip = sf_input("dip");
    out = sf_output("out");

    if (!sf_histint(dip,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(dip,"n2",&n2)) sf_error("No n2= in input");
    n12 = n1*n2;
    n3 = sf_leftsize(inp,2);

    if (!sf_getbool("verb",&verb)) verb=false;
    /* verbosity */
    if (!sf_getfloat("eps",&eps)) eps=0.01;
    /* regularization */
    
    if (!sf_getint("order",&order)) order=1;
    /* accuracy order */

    if (!sf_getint("rect",&n)) n=0;
    /* smoothing radius */

    predict_init (n1, n2, eps*eps, order, 1, false);

    p = sf_floatalloc2(n1,n2);
    u = sf_floatalloc2(n1,n2);

    next = sf_floatalloc2(n1,n2);
    trace = sf_floatalloc(n1);
    left = sf_floatalloc(n1);
    right =  sf_floatalloc(n1);

    for (i3=0; i3 < n3; i3++) {
	if (verb) fprintf(stderr,"cmp %d of %d\n",i3+1,n3);
	sf_floatread(u[0],n12,inp);
	sf_floatread(p[0],n12,dip);

	for (i=1; i < n; i++) {
	    r = 2*sinf(2*SF_PI*i/n);
	    r = 1/(r*r);
	
	    for (i2=0; i2 < n2; i2++) { /* loop over traces */	
		for  (i1=0; i1 < n1; i1++) {
		    trace[i1] = u[i2][i1];
		}

		/* prediction from the left */
		if (i2 > 0) {
		    for  (i1=0; i1 < n1; i1++) {
			left[i1] = u[i2-1][i1];
		    }
		    predict_step(false,true,left,p[i2-1]);
		} else {
		    for  (i1=0; i1 < n1; i1++) {
			left[i1] = trace[i1];
		    }
		}
		    
		/* prediction from the right */
		if (i2 < n2-1) {
		    for  (i1=0; i1 < n1; i1++) {
			right[i1] = u[i2+1][i1];
		    }
		    predict_step(false,false,right,p[i2]);
		} else {
		    for  (i1=0; i1 < n1; i1++) {
			right[i1] = trace[i1];
		    }
		}
		
		for  (i1=0; i1 < n1; i1++) {
		    next[i2][i1] = (1-2*r)*trace[i1]+r*(left[i1]+right[i1]);
		}
	    }

	    for (i2=0; i2 < n2; i2++) { 
		for  (i1=0; i1 < n1; i1++) {
		    u[i2][i1] = next[i2][i1];
		}
	    }
	}
	      	    
	sf_floatwrite(u[0],n12,out);
    }	    

    exit (0);
}

