/* 3-D plane-wave smoothing by box cascade. */
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
    bool verb, edge;
    int n1,n2,n3, n12, i1,i2,i3, order, i, n, nclip, ic, nc;
    float ***u, ***p, ***q, ***next, *left, *right, *top, *bottom, *trace, ***weight;
    float eps, r, pclip, a;
    sf_file inp, out, dip;

    sf_init(argc,argv);
    inp = sf_input("in");
    dip = sf_input("dip");
    out = sf_output("out");

    if (!sf_histint(inp,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(inp,"n2",&n2)) sf_error("No n2= in input");
    if (!sf_histint(inp,"n3",&n3)) sf_error("No n3= in input");
    n12 = n1*n2*n3;

    if (!sf_getbool("verb",&verb)) verb=false;
    /* verbosity */
    if (!sf_getfloat("eps",&eps)) eps=0.01;
    /* regularization */
    if (!sf_getbool("edge",&edge)) edge=false;
    /* preserve edges */
    
    if (!sf_getint("order",&order)) order=1;
    /* accuracy order */

    if (!sf_getint("rect",&n)) n=0;
    /* smoothing radius */

    if (!sf_getfloat("pclip",&pclip)) pclip=50.;
    /* percentage clip for the gradient */

    if (!sf_getint("cycle",&nc)) nc=1;
    /* number of cycles */
    
    nclip = (int) n12*pclip*0.01;
    if (nclip < 1) {
	nclip = 1;
    } else if (nclip > n12) {
	nclip = n12;
    }
    nclip--;

    predict_init (n1, n2, eps*eps, order, 1, false);

    p = sf_floatalloc3(n1,n2,n3);
    q = sf_floatalloc3(n1,n2,n3);
    u = sf_floatalloc3(n1,n2,n3);

    next = sf_floatalloc3(n1,n2,n3);
    trace =  sf_floatalloc(n1);
    left =   sf_floatalloc(n1);
    right =  sf_floatalloc(n1);
    top =    sf_floatalloc(n1);
    bottom = sf_floatalloc(n1);
    
    if (edge) {
	weight = sf_floatalloc3(n1,n2,n3);
	for (i3=0; i3 < n3; i3++) {
	    for (i2=0; i2 < n2; i2++) { 
		for  (i1=0; i1 < n1; i1++) {
		    weight[i3][i2][i1] = 1.0f;
		}
	    }
	}
    } else {
	weight = NULL;
    }

    sf_floatread(u[0][0],n12,inp);
    sf_floatread(p[0][0],n12,dip);
    sf_floatread(q[0][0],n12,dip);

    for (ic=0; ic < nc; ic++) {
	for (i=1; i < n; i++) {
	    sf_warning("cascade %d of %d (cycle %d)",i+1,n,ic+1);
	    
	    r = 2*sinf(2*SF_PI*i/n);
	    r = 0.5/(r*r);

	    for (i3=0; i3 < n3; i3++) {
		for (i2=0; i2 < n2; i2++) { /* loop over traces */	
		    for  (i1=0; i1 < n1; i1++) {
			trace[i1] = u[i3][i2][i1];
		    }
		    
		    /* prediction from the left */
		    if (i2 > 0) {
			for  (i1=0; i1 < n1; i1++) {
			    left[i1] = u[i3][i2-1][i1];
			}
			predict_step(false,true,left,p[i3][i2-1]);
			for  (i1=0; i1 < n1; i1++) {
			    left[i1] -= trace[i1];
			}
		    } else {
			for  (i1=0; i1 < n1; i1++) {
			    left[i1] = 0.0f;
			}
		    }
		    
		    /* prediction from the right */
		    if (i2 < n2-1) {
			for  (i1=0; i1 < n1; i1++) {
			    right[i1] = u[i3][i2+1][i1];
			}
			predict_step(false,false,right,p[i3][i2]);
			for  (i1=0; i1 < n1; i1++) {
			    right[i1] -= trace[i1];
			}
		    } else {
			for  (i1=0; i1 < n1; i1++) {
			    right[i1] = 0.0f;
			}
		    }

		    /* prediction from the top */
		    if (i3 > 0) {
			for  (i1=0; i1 < n1; i1++) {
			    top[i1] = u[i3-1][i2][i1];
			}
			predict_step(false,true,top,q[i3-1][i2]);
			for  (i1=0; i1 < n1; i1++) {
			    top[i1] -= trace[i1];
			}
		    } else {
			for  (i1=0; i1 < n1; i1++) {
			    top[i1] = 0.0f;
			}
		    }

		    /* prediction from the bottom */
		    if (i3 < n3-1) {
			for  (i1=0; i1 < n1; i1++) {
			    bottom[i1] = u[i3+1][i2][i1];
			}
			predict_step(false,false,bottom,q[i3][i2]);
			for  (i1=0; i1 < n1; i1++) {
			    bottom[i1] -= trace[i1];
			}
		    } else {
			for  (i1=0; i1 < n1; i1++) {
			    bottom[i1] = 0.0f;
			}
		    }
		 
		    if (edge) {
			for  (i1=0; i1 < n1; i1++) {
			    next[i3][i2][i1] = trace[i1]+ 
				r*weight[i3][i2][i1]*(left[i1]+right[i1]+top[i1]+bottom[i1]);
			    weight[i3][i2][i1] = 
				left[i1]*left[i1]+right[i1]*right[i1]+top[i1]*top[i1]+bottom[i1]*bottom[i1];
			}
		    } else {
			for  (i1=0; i1 < n1; i1++) {
			    next[i3][i2][i1] = trace[i1]+r*(left[i1]+right[i1]+top[i1]+bottom[i1]);
			}
		    }
		}
	    }

	    for (i3=0; i3 < n3; i3++) {	
		for (i2=0; i2 < n2; i2++) { 
		    for  (i1=0; i1 < n1; i1++) {
			u[i3][i2][i1] = next[i3][i2][i1];
		    }
		}
	    }

	    if (edge) {
		for (i3=0; i3 < n3; i3++) {	
		    for (i2=0; i2 < n2; i2++) { 
			for  (i1=0; i1 < n1; i1++) {
			    next[i3][i2][i1] = weight[i3][i2][i1];
			}
		    }
		}
		a = sf_quantile(nclip,n12,next[0][0]);
		for (i3=0; i3 < n3; i3++) {	
		    for (i2=0; i2 < n2; i2++) { 
			for  (i1=0; i1 < n1; i1++) {
			    weight[i3][i2][i1] = 1.0f/(1.0f+weight[i3][i2][i1]/a); 
			}
		    }
		}
	    }
	}
    }
	      	    
    sf_floatwrite(u[0][0],n12,out);

    exit (0);
}

