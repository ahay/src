/* Bilateral stacking. */
/*
  Copyright (C) 2009 University of Texas at Austin
  
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

int main (int argc, char *argv[])
{
    bool verb, bilat;
    int n1, n2, n3, n12, i1,i2,i3, foldw, foldn, iter, niter;
    float *w, *norm, *ref, *trace, *result, ax, bx, maxa, maxb;
    sf_file in, out, weight;

    sf_init(argc,argv);
    in  = sf_input("in");
    out = sf_output("out");
    weight = sf_output("weight");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    n12 = n1*n2;
    n3 = sf_leftsize(in,2);

    if (!sf_getbool("verb",&verb)) verb=false;
    /* verbosity */

    if (!sf_getint("niter",&niter)) niter=20;
    /* number of iterations */
    
    if (!sf_getbool("bilat",&bilat)) bilat=false;
    /* if y, bilateral smoothing */

    if (!sf_getfloat("ax",&ax)) sf_error("Need ax=");
    /* Gaussian weight for the range distance */
    
    if (!sf_getfloat("bx",&bx)) sf_error("Need bx=");
    /* Exponential weight for the domain distance */
    
    sf_unshiftdim(in, out, 2);	

    trace = sf_floatalloc(n12);
    ref = sf_floatalloc(n1);
    w = sf_floatalloc(n12);
    result = sf_floatalloc(n1);
    norm = sf_floatalloc(n1);

    for (i3=0; i3 < n3; i3++) {
	if (verb) sf_warning("slice %d of %d",i3+1,n3);

	sf_floatread(trace,n12,in);

	for (i1=0; i1 < n1; i1++) { /* initialize reference */
	    ref[i1] = 0.;
	}
	for (i1=0; i1 < n1; i1++) { /* reference stacking */
	    foldn = 0;
	    for (i2=0; i2 < n2; i2++) { 
		ref[i1] += trace[i2*n1+i1];
		if (0. != trace[i2*n1+i1]) foldn++;
	    }
	    ref[i1] = ref[i1]/(foldn+FLT_EPSILON);
	}
	   
	for(i1=0; i1 < n1; i1++) {
	    result[i1] = 0.;
	}	  
	for(iter=0; iter < niter; iter++) {

	    /* Scaling factor */
	    for(i1=0; i1 < n1; i1++) {
		for(i2=0; i2 < n2; i2++) {
		    w[i2*n1+i1] = trace[i2*n1+i1] - ref[i1];
		}
	    }
	    
	    maxa = 0.;
	    for(i1=0; i1 < n1; i1++) {
		if (maxa < fabsf(ref[i1])) maxa = fabsf(ref[i1]);
	    }
	    
	    maxb = 0.;
	    for(i1=0; i1 < n1; i1++) {
		for(i2=0; i2 < n2; i2++) {
		    if (maxb < fabsf(w[i2*n1+i1])) maxb = fabsf(w[i2*n1+i1]);
		}
	    }
	    
	    /* Define weights */
	    if (bilat) {
		for(i1=0; i1 < n1; i1++) {
		    for(i2=0; i2 < n2; i2++) {
			w[i2*n1+i1] = expf(ax*(ref[i1])*(ref[i1])
					   /(maxa*maxa+FLT_EPSILON)) * 
			    expf(-0.5*(trace[i2*n1+i1]-ref[i1])*(trace[i2*n1+i1]-ref[i1])
				 /(bx*bx*maxb*maxb+FLT_EPSILON));
		    }
		}
	    } else {
		for(i1=0; i1 < n1; i1++) {
		    for(i2=0; i2 < n2; i2++) {
			w[i2*n1+i1] = expf(-0.5*(trace[i2*n1+i1]-ref[i1])*(trace[i2*n1+i1]-ref[i1])
					   /(bx*bx*maxb*maxb+FLT_EPSILON));
		    }
		}
	    }
	    for(i1=0; i1 < n1; i1++) {
		result[i1] = 0.;
		norm[i1] = 0.;
	    }
	    
	    for(i1=0; i1 < n1; i1++) {
		foldw = 0;
		foldn = 0;
		for(i2=0; i2 < n2; i2++) {
		    result[i1] += trace[i2*n1+i1]*w[i2*n1+i1];
		    if (0 != trace[i2*n1+i1]*w[i2*n1+i1]) foldn++;
		}
		
		for(i2=0; i2 < n2; i2++) {
		    norm[i1] += w[i2*n1+i1];
		    if (0 != w[i2*n1+i1]) foldw++;
		}
		result[i1] = (result[i1]*foldw)/(norm[i1]*foldn+FLT_EPSILON);
	    }
	    
	    for(i1=0; i1 < n1; i1++) {
		ref[i1] = result[i1];
	    }	    
	}

	sf_floatwrite(w,n12,weight);	
	sf_floatwrite(result,n1,out);
    }	    

    exit (0);
}

/* 	$Id$	 */
