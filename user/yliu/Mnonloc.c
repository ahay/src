/* Non-local (Bilateral) smoothing. */
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
#include <stdio.h>
#include <math.h>

int main (int argc, char *argv[])
{
    int n1,n2, i1,i2, is, ns;
    float *trace, *trace2, ax, bx, t, norm;
    sf_file inp, out;

    /* initialize */
    sf_init(argc,argv);

    /* set input and output files */
    inp = sf_input("in");
    out = sf_output("out");

    /* get input dimensions */
    if (!sf_histint(inp,"n1",&n1)) 
	sf_error("No n1= in input");
    n2 = sf_leftsize(inp,1);
    
    /* get command-line parameters */
    if (!sf_getint("ns",&ns)) sf_error("Need ns=");
    /* spray radius */

    if (!sf_getfloat("ax",&ax)) sf_error("Need ax=");
    /* exponential weight for the range distance */

    if (!sf_getfloat("bx",&bx)) sf_error("Need bx=");
    /* exponential weight for the domain distance */
    trace = sf_floatalloc(n1);
    trace2 = sf_floatalloc(n1);

    /* loop over traces */
    for (i2=0; i2 < n2; i2++) {
	/* read input */
	sf_floatread(trace,n1,inp);

	/* loop over samples */
	for (i1=0; i1 < n1; i1++) {	    
	    t = 0.;
	    norm = 0.;

	    /* accumulate shifts */
	    for (is=-ns; is <= ns; is++) {
		if (i1+is >= 0 && i1+is < n1) {
		    t += trace[i1+is]*(1-is*is/(ax+FLT_EPSILON))*expf(-bx*(trace[i1+is]-trace[i1])*(trace[i1+is]-trace[i1]));
		    norm += (1-is*is/(ax+FLT_EPSILON))*expf(-bx*(trace[i1+is]-trace[i1])*(trace[i1+is]-trace[i1]));
/*		    t += trace[i1+is]*expf(-0.5*is*is/(ax*ax+FLT_EPSILON))*expf(-0.5*(trace[i1+is]-trace[i1])*(trace[i1+is]-trace[i1])/(bx*bx+FLT_EPSILON));
		    norm += expf(-0.5*is*is/(ax*ax+FLT_EPSILON))*expf(-0.5*(trace[i1+is]-trace[i1])*(trace[i1+is]-trace[i1])/(bx*bx+FLT_EPSILON));
*/		} 
	    }
	    /* Nomalize */
	    trace2[i1] = t / (norm+FLT_EPSILON);
	}

	/* write output */
	sf_floatwrite(trace2,n1,out);
    }

    /* clean up */
    sf_fileclose(inp);
    exit (0);
}

/* 	$Id$	 */
