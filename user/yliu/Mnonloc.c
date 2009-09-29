/* Non-local (Bilateral) smoothing. 
Tomasi, C., and R. Manduchi, 1998, Bilateral filtering 
for gray and color images: Proceedings of the 1998 
IEEE International Conference on Computer Vision, 
IEEE, 836-846*/
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
    int n1,n2, i1,i2, is, ns, irep, nrep;
    float *trace, *output, ax, bx, t, norm;
    sf_file inp, out;
    bool gauss;

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

    if (!sf_getfloat("bx",&bx)) sf_error("Need bx=");
    /* exponential weight for the domain distance (different if gaussian)*/

    if (!sf_getbool("gauss",&gauss)) gauss=false;
    /* if y, Gaussian weight, whereas Triangle weight */

    if (!sf_getint("repeat",&nrep)) nrep=1;
    /* repeat filtering several times */

    if (gauss) {
	if (!sf_getfloat("ax",&ax)) sf_error("Need ax=");
	/* Gaussian weight for the range distance */
    }

    trace = sf_floatalloc(n1);
    output = sf_floatalloc(n1);

    /* loop over traces */
    for (i2=0; i2 < n2; i2++) {
	/* read input */
	sf_floatread(trace,n1,inp);
	for (i1=0; i1 < n1; i1++) {
	    output[i1] = trace[i1];
	}

	for (irep=0; irep < nrep; irep++) {
	    /* loop over samples */
	    for (i1=0; i1 < n1; i1++) {	    
		t = 0.;
		norm = 0.;
		
		/* accumulate shifts */
		for (is=-ns; is <= ns; is++) {
		    if (i1+is >= 0 && i1+is < n1) {
			if (gauss) {
			    t += output[i1+is]*expf(-0.5*is*is/(ax*ax+FLT_EPSILON))*expf(-0.5*(output[i1+is]-output[i1])*(output[i1+is]-output[i1])/(bx*bx+FLT_EPSILON));
			    norm += expf(-0.5*is*is/(ax*ax+FLT_EPSILON))*expf(-0.5*(output[i1+is]-output[i1])*(output[i1+is]-output[i1])/(bx*bx+FLT_EPSILON));
			} else {
			    t += output[i1+is]*(1.-fabsf(1.*is)/(ns+FLT_EPSILON))*expf(-bx*(output[i1+is]-output[i1])*(output[i1+is]-output[i1]));
			    norm += (1.-fabsf(1.*is)/(ns+FLT_EPSILON))*expf(-bx*(output[i1+is]-output[i1])*(output[i1+is]-output[i1]));
			}
		    } 
		}
		/* Nomalize */
		output[i1] = t / (norm+FLT_EPSILON);
	    }
	}

	/* write output */
	sf_floatwrite(output,n1,out);
    }

    /* clean up */

    exit (0);
}

/* 	$Id$	 */
